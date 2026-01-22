import express from "express";
import cors from "cors";
import crypto from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";
import { exec } from "node:child_process";
import { fileURLToPath } from "node:url";
import {
  query,
  isSDKAssistantMessage,
  isSDKPartialAssistantMessage,
  isSDKResultMessage,
  type CanUseTool,
} from "@qwen-code/sdk";

// 定义 __dirname（ESM 模块中需要）
const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const WORKDIR = "D:\\threo333\\Projects\\AAA_delta_GPA\\codes";
const QWEN_EXE = "node:C:\\Users\\threo\\AppData\\Roaming\\npm\\node_modules\\@qwen-code\\qwen-code\\cli.js";
const QWEN_MD_PATH = path.join(WORKDIR, "QWEN.md");

// === 会话与历史（内存版） ===
// 生产环境可换成 sqlite/文件/redis；本机自用内存足够
type ChatMsg = { role: "user" | "assistant"; content: string; ts: number };
const histories = new Map<string, ChatMsg[]>();

function getHistory(sessionId: string) {
  let h = histories.get(sessionId);
  if (!h) {
    h = [];
    histories.set(sessionId, h);
  }
  return h;
}

function pushHistory(sessionId: string, msg: ChatMsg) {
  const h = getHistory(sessionId);
  h.push(msg);
  // 简单限长：最多保留 40 条消息（20 轮）
  const MAX = 40;
  if (h.length > MAX) h.splice(0, h.length - MAX);
}

function formatHistoryForPrompt(h: ChatMsg[]) {
  // 用清晰的对话格式（对模型更稳）
  return h
    .map((m) => (m.role === "user" ? `用户：${m.content}` : `助手：${m.content}`))
    .join("\n\n");
}

// === QWEN.md 缓存 ===
let qwenMdText: string | null = null;
async function getQwenMdText() {
  if (qwenMdText != null) return qwenMdText;
  qwenMdText = await fs.readFile(QWEN_MD_PATH, "utf-8");
  return qwenMdText;
}

// === SSE 管理 ===
type SseClient = {
  res: express.Response;
  send: (event: string, data: unknown) => void;
};

const sseClients = new Map<string, SseClient>();

function makeSseSender(res: express.Response) {
  return (event: string, data: unknown) => {
    try {
      res.write(`event: ${event}\n`);
      res.write(`data: ${JSON.stringify(data)}\n\n`);
    } catch (e) {
      // 连接可能已断开，忽略写入错误
    }
  };
}

// 动态获取并发送SSE消息的辅助函数（每次发送前重新获取最新的client）
function sendToSession(sessionId: string, event: string, data: unknown): boolean {
  const client = sseClients.get(sessionId);
  if (client) {
    client.send(event, data);
    return true;
  }
  return false;
}

// 记录已推送过工作区事件的会话，避免重复推送
const sessionWorkspaceNotified = new Set<string>();

// 检测工作区创建（通过工具调用参数）
// 当检测到创建 works/YYYYMMDD_XXX_... 格式的文件夹时，自动设置工作区
function detectWorkspaceCreation(sessionId: string, toolName: string, input: any) {
  // 如果该会话已经推送过工作区事件，跳过
  if (sessionWorkspaceNotified.has(sessionId)) return;
  
  try {
    // 解析 input
    const inputObj = typeof input === 'string' ? JSON.parse(input) : input;
    if (!inputObj) return;
    
    let targetPath: string | null = null;
    
    // 方式1：检测 shell 命令（如 run_shell_command, run_in_terminal 等）
    if (inputObj.command && typeof inputObj.command === 'string') {
      const cmd = inputObj.command;
      // 匹配 mkdir 命令：mkdir works\YYYYMMDD_XXX_name 或 mkdir -p works/...
      const mkdirPatterns = [
        /mkdir\s+(?:-p\s+)?["']?([^"'\s]+works[\\\/]\d{8}_\d{3}_[\w\-]+)/i,
        /mkdir\s+(?:-p\s+)?["']?.*?(works[\\\/]\d{8}_\d{3}_[\w\-]+)/i,
        /New-Item\s+.*?["']?([^"'\s]*works[\\\/]\d{8}_\d{3}_[\w\-]+)/i,
      ];
      
      for (const pattern of mkdirPatterns) {
        const match = cmd.match(pattern);
        if (match) {
          targetPath = match[1];
          break;
        }
      }
    }
    
    // 方式2：检查常见的文件/目录创建工具的路径字段
    if (!targetPath) {
      const pathFields = ['path', 'filePath', 'dirPath', 'directory', 'folder'];
      for (const field of pathFields) {
        if (inputObj[field] && typeof inputObj[field] === 'string') {
          targetPath = inputObj[field];
          break;
        }
      }
    }
    
    if (!targetPath) return;
    
    // 匹配 works/YYYYMMDD_XXX_name 格式
    const workspacePattern = /works[\\\/](\d{8}_\d{3}_[\w\-]+)/;
    const match = targetPath.match(workspacePattern);
    
    if (match) {
      const workspaceName = match[1]; // 例如 "20260122_001_newtask"
      const workspacePath = `works/${workspaceName}`;
      
      console.log("[Workspace] 检测到工作区创建:", { toolName, workspacePath, command: inputObj.command });
      
      // 标记已通知，避免重复推送
      sessionWorkspaceNotified.add(sessionId);
      
      // 自动设置该会话的工作区
      sessionWorkspaces.set(sessionId, path.join(WORKDIR, workspacePath));
      
      // 推送工作区事件给前端
      sendToSession(sessionId, "workspace", { 
        path: workspacePath,
        absolutePath: path.join(WORKDIR, workspacePath)
      });
    }
  } catch (e) {
    // 忽略解析错误
  }
}

// === 工具审批 pending 表 ===
type PendingApproval = {
  resolve: (r: { behavior: "allow"; updatedInput: any } | { behavior: "deny"; message?: string }) => void;
  timer: NodeJS.Timeout;
  sessionId: string;
};

const pendingApprovals = new Map<string, PendingApproval>();

// === 会话中止控制器 ===
const sessionAbortControllers = new Map<string, AbortController>();

// === App ===
const app = express();
app.use(cors());
app.use(express.json({ limit: "2mb" }));

// 提供静态文件（前端）
const distPath = path.join(__dirname, "dist");
app.use(express.static(distPath));

// 请求日志
app.use((req, res, next) => {
  console.log(`[HTTP] ${req.method} ${req.path}`);
  next();
});

app.get("/api/sse", (req, res) => {
  const sessionId = String(req.query.sessionId ?? "").trim();
  if (!sessionId) return res.status(400).json({ error: "missing sessionId" });

  console.log("[SSE] connect", { sessionId, ip: req.ip });

  res.status(200);
  res.setHeader("Content-Type", "text/event-stream; charset=utf-8");
  res.setHeader("Cache-Control", "no-cache, no-transform");
  res.setHeader("Connection", "keep-alive");
  // 立即发送 HTTP 头，防止某些浏览器/代理等待
  if (res.flushHeaders) res.flushHeaders();

  const send = makeSseSender(res);
  sseClients.set(sessionId, { res, send });

  send("meta", { ok: true, sessionId });

  req.on("close", () => {
    console.log("[SSE] close", { sessionId });
    sseClients.delete(sessionId);
  });

  req.on("error", (err) => {
    console.log("[SSE] req error", { sessionId, err: err?.message });
    sseClients.delete(sessionId);
  });

  res.on("error", (err) => {
    console.log("[SSE] res error", { sessionId, err: err?.message });
    sseClients.delete(sessionId);
  });
});


// 清空某个会话历史（相当于“新对话”）
app.post("/api/session/reset", (req, res) => {
  console.log("[CHAT] incoming", { sessionId: req.body?.sessionId });
  const { sessionId } = req.body ?? {};
  if (!sessionId) return res.status(400).json({ error: "missing sessionId" });
  histories.set(String(sessionId), []);  // 清除工作区状态
  sessionWorkspaces.delete(sessionId);
  sessionWorkspaceNotified.delete(sessionId);  res.json({ ok: true });
});

// 前端提交工具审批结果
app.post("/api/tool/decision", (req, res) => {
  res.json({ ok: true, message: "工具自动执行，无需审批" });
});

// 中止当前生成
app.post("/api/chat/interrupt", (req, res) => {
  const sessionId = String(req.body?.sessionId ?? "").trim();
  if (!sessionId) return res.status(400).json({ error: "missing sessionId" });

  const controller = sessionAbortControllers.get(sessionId);
  if (controller) {
    console.log("[CHAT] 用户中止生成", { sessionId });
    controller.abort();
    sessionAbortControllers.delete(sessionId);
    
    // 通知前端中止完成
    sendToSession(sessionId, "interrupted", { sessionId });
    res.json({ ok: true, message: "已中止" });
  } else {
    res.json({ ok: true, message: "没有正在进行的生成" });
  }
});

// 打开本地文件（用 VS Code 或系统默认程序）
app.post("/api/open-file", async (req, res) => {
  let filePath = String(req.body?.filePath ?? "").trim();
  let contextDir = String(req.body?.contextDir ?? "").trim(); // 上下文目录（当前工作文件夹）
  if (!filePath) return res.status(400).json({ error: "missing filePath" });

  console.log("[OPEN-FILE] 收到请求:", { filePath, contextDir });

  // 处理路径解析
  let normalizedPath: string;
  
  if (path.isAbsolute(filePath)) {
    // 已经是绝对路径
    normalizedPath = path.normalize(filePath);
    console.log("[OPEN-FILE] 使用绝对路径:", normalizedPath);
  } else if (contextDir && !filePath.includes('/') && !filePath.includes('\\')) {
    // 纯文件名 + 有上下文目录：优先在上下文目录中查找
    const contextPath = path.isAbsolute(contextDir) 
      ? path.normalize(path.join(contextDir, filePath))
      : path.normalize(path.join(WORKDIR, contextDir, filePath));
    normalizedPath = contextPath;
    console.log("[OPEN-FILE] 使用上下文目录:", normalizedPath);
  } else {
    // 相对路径，基于 WORKDIR 解析
    normalizedPath = path.normalize(path.join(WORKDIR, filePath));
    console.log("[OPEN-FILE] 使用 WORKDIR 相对路径:", normalizedPath);
  }

  // 安全检查：只允许打开工作目录下的文件
  if (!normalizedPath.toLowerCase().startsWith(WORKDIR.toLowerCase())) {
    console.log("[OPEN-FILE] 拒绝打开工作目录外的文件:", normalizedPath);
    return res.status(403).json({ error: "只能打开工作目录内的文件" });
  }

  try {
    // 检查文件是否存在
    await fs.access(normalizedPath);
    
    // 使用 VS Code 打开文件
    const command = `code "${normalizedPath}"`;
    console.log("[OPEN-FILE] 执行命令:", command);
    
    exec(command, (error) => {
      if (error) {
        console.error("[OPEN-FILE] 执行失败:", error.message);
        // 如果 VS Code 失败，尝试用系统默认程序
        exec(`start "" "${normalizedPath}"`, (err2) => {
          if (err2) {
            console.error("[OPEN-FILE] 系统打开也失败:", err2.message);
          }
        });
      }
    });
    
    res.json({ ok: true, path: normalizedPath });
  } catch (err: any) {
    console.error("[OPEN-FILE] 文件不存在:", normalizedPath);
    res.status(404).json({ error: "文件不存在", path: normalizedPath });
  }
});

// 发起一轮对话：前端 POST /api/chat {sessionId, prompt}
// 实际输出走 SSE（delta/assistant/done）
app.post("/api/chat", async (req, res) => {
  console.log("[CHAT-DEBUG-1] POST /api/chat 收到请求");
  
  const sessionId = String(req.body?.sessionId ?? "").trim();
  const userPrompt = String(req.body?.prompt ?? "").trim();
  
  console.log(`[CHAT-DEBUG-2] sessionId=${sessionId}, prompt=${userPrompt}`);
  
  if (!sessionId || !userPrompt) {
    console.log("[CHAT-DEBUG-3] 缺少参数");
    return res.status(400).json({ error: "missing sessionId or prompt" });
  }

  console.log("[CHAT] incoming", { sessionId, promptLen: userPrompt.length });

  // 立刻返回，避免前端卡住；流式输出通过 SSE 推送
  res.json({ ok: true });

  // 检查SSE连接是否存在（但不持有引用，每次发送时动态获取）
  if (!sseClients.has(sessionId)) {
    console.warn("[CHAT] SSE client not connected", { sessionId });
    return;
  }
  
  console.log("[CHAT-DEBUG-4] 开始处理，SSE client 已连接");

  // 记录用户消息到历史
  pushHistory(sessionId, { role: "user", content: userPrompt, ts: Date.now() });

  // === 构造“始终带着 QWEN.md + 多轮历史”的 prompt ===
  let qwenMd = "";
  try {
    qwenMd = await getQwenMdText();
  } catch (e: any) {    console.error("[CHAT] 读取 QWEN.md 失败：", e);    sendToSession(sessionId, "error", { message: `读取 QWEN.md 失败：${e?.message ?? String(e)}` });
    return;
  }

  // 防止 QWEN.md 或历史太长（你可按需要调大/调小）
  const MAX_QWEN_MD_CHARS = 20000;
  const qwenMdSafe =
    qwenMd.length > MAX_QWEN_MD_CHARS ? qwenMd.slice(0, MAX_QWEN_MD_CHARS) + "\n\n[...QWEN.md 已截断...]" : qwenMd;

  const history = getHistory(sessionId);
  console.log("[CHAT-DEBUG-HISTORY] 当前历史条数:", history.length, "内容摘要:", history.map(h => ({ role: h.role, len: h.content.length, preview: h.content.slice(0, 50) })));
  const formattedHistory = formatHistoryForPrompt(history.slice(0, -1)); // 不重复包含本轮用户消息（我们单独放在末尾）
  console.log("[CHAT-DEBUG-HISTORY] 格式化历史长度:", formattedHistory.length);

  const wrappedPrompt =
    `你是本项目的代码助手。QWEN.md 是最高优先级的项目说明与约束，必须严格遵守。\n\n` +
    `=== QWEN.md BEGIN ===\n${qwenMdSafe}\n=== QWEN.md END ===\n\n` +
    (formattedHistory ? `=== 对话历史 BEGIN ===\n${formattedHistory}\n=== 对话历史 END ===\n\n` : "") +
    `用户：${userPrompt}\n\n` +
    `助手：`;

  sendToSession(sessionId, "start", { sessionId });
  console.log("[CHAT] 开始生成", { sessionId });

  // 使用 AbortController 实现超时控制和手动中止
  const abortController = new AbortController();
  sessionAbortControllers.set(sessionId, abortController);
  
  const timeout = setTimeout(() => {
    console.warn("[CHAT] 生成超时（10分钟），中止", { sessionId });
    abortController.abort();
  }, 600_000); // 10 分钟超时

  try {
    const result = query({
      prompt: wrappedPrompt,
      options: {
        cwd: WORKDIR,
        pathToQwenExecutable: QWEN_EXE,
        permissionMode: "yolo", // 让所有工具自动执行（包括文件修改和 Python 命令）
        includePartialMessages: true,
        abortController,
      },
    });

    let assistantFull = "";
    let messageCount = 0;
    
    // 记录本轮对话中的工具调用（用于生成历史摘要）
    const toolCallsForHistory: { name: string; input: string; result: string }[] = [];

    // 记录 streaming 中出现的 tool_use/tool_result（用于前端展示）
    const blockStates = new Map<number, { type: string; toolUseId?: string; name?: string; input?: string; result?: string; isError?: boolean }>();

    for await (const msg of result) {
      try {
        const msgObj = msg as any;
        
        if (msgObj?.type === "stream_event") {
          const evt = msgObj.event;

          // 标记 content_block 类型，便于后续增量拼装
          if (evt?.type === "content_block_start") {
            const cb = evt.content_block as any;
            blockStates.set(evt.index, {
              type: cb?.type,
              toolUseId: cb?.type === "tool_use" ? cb.id : cb?.tool_use_id,
              name: cb?.name,
              input: typeof cb?.input === "string" ? cb.input : undefined,
              result: typeof cb?.content === "string" ? cb.content : undefined,
              isError: cb?.is_error ?? false,
            });

            // 立即把开始事件推给前端
            if (cb?.type === "tool_use") {
              sendToSession(sessionId, "tool_event", {
                kind: "use",
                toolUseId: cb.id,
                name: cb.name,
                input: cb.input,
              });
              
              // 检测是否是创建工作区文件夹的操作
              detectWorkspaceCreation(sessionId, cb.name, cb.input);
            }
            if (cb?.type === "tool_result") {
              sendToSession(sessionId, "tool_event", {
                kind: "result",
                toolUseId: cb.tool_use_id,
                isError: cb.is_error ?? false,
                result: cb.content,
              });
            }
          }

          // 增量内容（text / json）
          if (evt?.type === "content_block_delta") {
            const state = blockStates.get(evt.index);
            if (state?.type === "tool_use" && evt.delta?.type === "input_json_delta") {
              const partial = evt.delta?.partial_json ?? "";
              state.input = (state.input ?? "") + partial;
              sendToSession(sessionId, "tool_event", {
                kind: "use",
                toolUseId: state.toolUseId,
                name: state.name,
                input: state.input,
                partial: true,
              });
              
              // 尝试检测工作区创建（input 可能已经包含完整的 command）
              detectWorkspaceCreation(sessionId, state.name || "", state.input);
            }
            if (state?.type === "tool_result" && evt.delta?.type === "text_delta") {
              const partial = evt.delta?.text ?? "";
              state.result = (state.result ?? "") + partial;
              sendToSession(sessionId, "tool_event", {
                kind: "result",
                toolUseId: state.toolUseId,
                result: state.result,
                isError: state.isError ?? false,
                partial: true,
              });
            }
          }

          // 处理文本流式输出
          const deltaText = evt?.delta?.text;
          if (typeof deltaText === "string" && deltaText.length > 0) {
            assistantFull += deltaText;
            sendToSession(sessionId, "delta", { text: deltaText });
            messageCount++;
          }
        } else if (msgObj?.type === "assistant") {
          console.log("[CHAT-ASSISTANT] received final message");
        } else if (msgObj?.type === "result") {
          let resultText = msgObj?.result;
          if (typeof resultText === "string" && resultText.length > assistantFull.length) {
            console.log("[CHAT-RESULT] using result as fallback, length:", resultText.length);
            assistantFull = resultText;
            sendToSession(sessionId, "delta", { text: resultText });
            messageCount++;
          }
        }
        // 其它消息类型先忽略
      } catch (msgErr: any) {
        console.log("[CHAT] 消息处理错误:", msgErr?.message);
      }
    }

    // 从 blockStates 中整理工具调用信息，生成历史摘要
    const toolSummaries: string[] = [];
    for (const [, state] of blockStates) {
      if (state.type === "tool_use" && state.name) {
        // 简化 input（截断过长的内容）
        let inputSummary = "";
        try {
          const inputObj = typeof state.input === 'string' ? JSON.parse(state.input) : state.input;
          inputSummary = JSON.stringify(inputObj);
          if (inputSummary.length > 200) inputSummary = inputSummary.slice(0, 200) + "...";
        } catch {
          inputSummary = state.input?.slice(0, 200) || "";
        }
        toolSummaries.push(`[工具调用: ${state.name}] ${inputSummary}`);
      }
      if (state.type === "tool_result" && state.result) {
        // 简化 result（截断过长的内容）
        const resultSummary = state.result.length > 300 ? state.result.slice(0, 300) + "..." : state.result;
        toolSummaries.push(`[工具结果${state.isError ? '(错误)' : ''}] ${resultSummary}`);
      }
    }
    
    // 将工具调用摘要附加到助手回复
    let contentForHistory = assistantFull.trim();
    if (toolSummaries.length > 0) {
      contentForHistory += "\n\n--- 本轮工具调用记录 ---\n" + toolSummaries.join("\n");
    }

    // 把本轮助手最终内容写入历史（包含工具调用摘要）
    if (contentForHistory) {
      console.log("[CHAT-DEBUG-SAVE] 保存助手回复，长度:", contentForHistory.length, "预览:", contentForHistory.slice(0, 100), "工具调用数:", toolSummaries.length);
      pushHistory(sessionId, { role: "assistant", content: contentForHistory, ts: Date.now() });
      sendToSession(sessionId, "done", { ok: true });
      console.log("[CHAT] 生成完成", { sessionId, responseLen: contentForHistory.length, messageCount, toolCalls: toolSummaries.length });
    } else {
      console.warn("[CHAT] 空响应", { sessionId });
      sendToSession(sessionId, "error", { message: "生成的回复为空，请检查 Qwen Code 配置" });
    }
  } catch (err: any) {
    console.error("[CHAT] 生成错误:", {
      sessionId,
      error: err?.message ?? String(err),
      stack: err?.stack,
    });
    
    // 向前端发送详细的错误信息
    const errorMsg = err?.message ?? String(err);
    const userFriendlyMsg = errorMsg.includes("ENOENT")
      ? `找不到 Qwen 可执行文件: ${QWEN_EXE}`
      : errorMsg.includes("ETIMEDOUT")
      ? "请求超时，Qwen Code 无响应"
      : errorMsg.includes("aborted")
      ? "生成超时（10分钟）"
      : errorMsg;
    
    sendToSession(sessionId, "error", { 
      message: `[Qwen 错误] ${userFriendlyMsg}`,
      debug: process.env.DEBUG ? errorMsg : undefined 
    });
  } finally {
    clearTimeout(timeout);
    sessionAbortControllers.delete(sessionId);
  }
});

// === 工作区管理 ===
// 存储每个会话的当前工作区路径
const sessionWorkspaces = new Map<string, string>();

// 设置当前工作区
app.post("/api/workspace/set", (req, res) => {
  const sessionId = String(req.body?.sessionId ?? "").trim();
  const workspacePath = String(req.body?.workspacePath ?? "").trim();
  
  if (!sessionId) return res.status(400).json({ error: "missing sessionId" });
  if (!workspacePath) return res.status(400).json({ error: "missing workspacePath" });
  
  // 构建完整路径
  const fullPath = path.isAbsolute(workspacePath) 
    ? workspacePath 
    : path.join(WORKDIR, workspacePath);
  
  // 安全检查
  if (!fullPath.toLowerCase().startsWith(WORKDIR.toLowerCase())) {
    return res.status(403).json({ error: "工作区必须在项目目录内" });
  }
  
  sessionWorkspaces.set(sessionId, fullPath);
  console.log("[WORKSPACE] 设置工作区:", { sessionId, path: fullPath });
  res.json({ ok: true, path: fullPath });
});

// 获取当前工作区
app.get("/api/workspace/get", (req, res) => {
  const sessionId = String(req.query.sessionId ?? "").trim();
  if (!sessionId) return res.status(400).json({ error: "missing sessionId" });
  
  const workspacePath = sessionWorkspaces.get(sessionId);
  res.json({ ok: true, path: workspacePath ?? null });
});

// 获取工作区文件列表（递归）
app.get("/api/workspace/files", async (req, res) => {
  const sessionId = String(req.query.sessionId ?? "").trim();
  if (!sessionId) return res.status(400).json({ error: "missing sessionId" });
  
  const workspacePath = sessionWorkspaces.get(sessionId);
  if (!workspacePath) {
    return res.json({ ok: true, files: [], message: "未设置工作区" });
  }
  
  try {
    await fs.access(workspacePath);
  } catch {
    return res.status(404).json({ error: "工作区目录不存在" });
  }
  
  interface FileNode {
    name: string;
    path: string;
    relativePath: string;
    isDirectory: boolean;
    children?: FileNode[];
  }
  
  async function scanDir(dirPath: string, relativePath: string = ""): Promise<FileNode[]> {
    const entries = await fs.readdir(dirPath, { withFileTypes: true });
    const result: FileNode[] = [];
    
    // 排序：文件夹在前，然后按名称排序
    entries.sort((a, b) => {
      if (a.isDirectory() && !b.isDirectory()) return -1;
      if (!a.isDirectory() && b.isDirectory()) return 1;
      return a.name.localeCompare(b.name);
    });
    
    for (const entry of entries) {
      // 跳过隐藏文件和常见忽略目录
      if (entry.name.startsWith('.') || 
          entry.name === 'node_modules' || 
          entry.name === '__pycache__' ||
          entry.name === '.ipynb_checkpoints') {
        continue;
      }
      
      const fullPath = path.join(dirPath, entry.name);
      const relPath = relativePath ? `${relativePath}/${entry.name}` : entry.name;
      
      if (entry.isDirectory()) {
        const children = await scanDir(fullPath, relPath);
        result.push({
          name: entry.name,
          path: fullPath,
          relativePath: relPath,
          isDirectory: true,
          children
        });
      } else {
        result.push({
          name: entry.name,
          path: fullPath,
          relativePath: relPath,
          isDirectory: false
        });
      }
    }
    
    return result;
  }
  
  try {
    const files = await scanDir(workspacePath);
    res.json({ ok: true, files, workspacePath });
  } catch (err: any) {
    console.error("[WORKSPACE] 扫描文件失败:", err);
    res.status(500).json({ error: "扫描文件失败: " + err.message });
  }
});

// 上传文件到工作区
import multer from "multer";

const upload = multer({ 
  storage: multer.memoryStorage(),
  limits: { fileSize: 100 * 1024 * 1024 } // 100MB 限制
});

app.post("/api/workspace/upload", upload.array("files"), async (req, res) => {
  const sessionId = String(req.body?.sessionId ?? "").trim();
  if (!sessionId) return res.status(400).json({ error: "missing sessionId" });
  
  const workspacePath = sessionWorkspaces.get(sessionId);
  if (!workspacePath) {
    return res.status(400).json({ error: "未设置工作区" });
  }
  
  const files = req.files as Express.Multer.File[];
  if (!files || files.length === 0) {
    return res.status(400).json({ error: "没有上传文件" });
  }
  
  const results: { name: string; success: boolean; error?: string }[] = [];
  
  for (const file of files) {
    const targetPath = path.join(workspacePath, file.originalname);
    try {
      await fs.writeFile(targetPath, file.buffer);
      results.push({ name: file.originalname, success: true });
      console.log("[WORKSPACE] 文件上传成功:", targetPath);
    } catch (err: any) {
      results.push({ name: file.originalname, success: false, error: err.message });
      console.error("[WORKSPACE] 文件上传失败:", err);
    }
  }
  
  res.json({ ok: true, results });
});

// 用系统默认程序打开工作区中的文件
app.post("/api/workspace/open", async (req, res) => {
  const sessionId = String(req.body?.sessionId ?? "").trim();
  const relativePath = String(req.body?.relativePath ?? "").trim();
  
  if (!sessionId) return res.status(400).json({ error: "missing sessionId" });
  if (!relativePath) return res.status(400).json({ error: "missing relativePath" });
  
  const workspacePath = sessionWorkspaces.get(sessionId);
  if (!workspacePath) {
    return res.status(400).json({ error: "未设置工作区" });
  }
  
  const fullPath = path.join(workspacePath, relativePath);
  
  // 安全检查
  if (!fullPath.toLowerCase().startsWith(workspacePath.toLowerCase())) {
    return res.status(403).json({ error: "不能打开工作区外的文件" });
  }
  
  try {
    await fs.access(fullPath);
    
    // 使用 VS Code 打开
    exec(`code "${fullPath}"`, (error) => {
      if (error) {
        // 如果 VS Code 失败，用系统默认程序
        exec(`start "" "${fullPath}"`);
      }
    });
    
    res.json({ ok: true, path: fullPath });
  } catch {
    res.status(404).json({ error: "文件不存在" });
  }
});

// 只监听本机，避免局域网其他设备访问
app.listen(915, "127.0.0.1", () => {
  console.log("\n====================================");
  console.log("✓ Server running: http://127.0.0.1:915");
  console.log("WORKDIR:", WORKDIR);
  console.log("QWEN_EXE:", QWEN_EXE);
  console.log("====================================\n");
});
