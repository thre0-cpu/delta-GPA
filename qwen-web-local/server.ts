import express from "express";
import cors from "cors";
import crypto from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";
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
    res.write(`event: ${event}\n`);
    res.write(`data: ${JSON.stringify(data)}\n\n`);
  };
}

// === 工具审批 pending 表 ===
type PendingApproval = {
  resolve: (r: { behavior: "allow"; updatedInput: any } | { behavior: "deny"; message?: string }) => void;
  timer: NodeJS.Timeout;
  sessionId: string;
};

const pendingApprovals = new Map<string, PendingApproval>();

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
  histories.set(String(sessionId), []);
  res.json({ ok: true });
});

// 前端提交工具审批结果
app.post("/api/tool/decision", (req, res) => {
  res.json({ ok: true, message: "工具自动执行，无需审批" });
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

  const client = sseClients.get(sessionId);
  if (!client) {
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
  } catch (e: any) {    console.error("[CHAT] 读取 QWEN.md 失败：", e);    client.send("error", { message: `读取 QWEN.md 失败：${e?.message ?? String(e)}` });
    return;
  }

  // 防止 QWEN.md 或历史太长（你可按需要调大/调小）
  const MAX_QWEN_MD_CHARS = 20000;
  const qwenMdSafe =
    qwenMd.length > MAX_QWEN_MD_CHARS ? qwenMd.slice(0, MAX_QWEN_MD_CHARS) + "\n\n[...QWEN.md 已截断...]" : qwenMd;

  const history = getHistory(sessionId);
  const formattedHistory = formatHistoryForPrompt(history.slice(0, -1)); // 不重复包含本轮用户消息（我们单独放在末尾）

  const wrappedPrompt =
    `你是本项目的代码助手。QWEN.md 是最高优先级的项目说明与约束，必须严格遵守。\n\n` +
    `=== QWEN.md BEGIN ===\n${qwenMdSafe}\n=== QWEN.md END ===\n\n` +
    (formattedHistory ? `=== 对话历史 BEGIN ===\n${formattedHistory}\n=== 对话历史 END ===\n\n` : "") +
    `用户：${userPrompt}\n\n` +
    `助手：`;

  client.send("start", { sessionId });
  console.log("[CHAT] 开始生成", { sessionId });

  // 使用 AbortController 实现超时控制
  const abortController = new AbortController();
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
              client.send("tool_event", {
                kind: "use",
                toolUseId: cb.id,
                name: cb.name,
                input: cb.input,
              });
            }
            if (cb?.type === "tool_result") {
              client.send("tool_event", {
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
              client.send("tool_event", {
                kind: "use",
                toolUseId: state.toolUseId,
                name: state.name,
                input: state.input,
                partial: true,
              });
            }
            if (state?.type === "tool_result" && evt.delta?.type === "text_delta") {
              const partial = evt.delta?.text ?? "";
              state.result = (state.result ?? "") + partial;
              client.send("tool_event", {
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
            client.send("delta", { text: deltaText });
            messageCount++;
          }
        } else if (msgObj?.type === "assistant") {
          console.log("[CHAT-ASSISTANT] received final message");
        } else if (msgObj?.type === "result") {
          let resultText = msgObj?.result;
          if (typeof resultText === "string" && resultText.length > assistantFull.length) {
            console.log("[CHAT-RESULT] using result as fallback, length:", resultText.length);
            assistantFull = resultText;
            client.send("delta", { text: resultText });
            messageCount++;
          }
        }
        // 其它消息类型先忽略
      } catch (msgErr: any) {
        console.log("[CHAT] 消息处理错误:", msgErr?.message);
      }
    }

    // 把本轮助手最终内容写入历史（用 assistantFull，保证和网页上看到的一致）
    if (assistantFull.trim()) {
      pushHistory(sessionId, { role: "assistant", content: assistantFull.trim(), ts: Date.now() });
      client.send("done", { ok: true });
      console.log("[CHAT] 生成完成", { sessionId, responseLen: assistantFull.length, messageCount });
    } else {
      console.warn("[CHAT] 空响应", { sessionId });
      client.send("error", { message: "生成的回复为空，请检查 Qwen Code 配置" });
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
    
    client.send("error", { 
      message: `[Qwen 错误] ${userFriendlyMsg}`,
      debug: process.env.DEBUG ? errorMsg : undefined 
    });
  } finally {
    clearTimeout(timeout);
  }
});

// 只监听本机，避免局域网其他设备访问
app.listen(3000, "127.0.0.1", () => {
  console.log("\n====================================");
  console.log("✓ Server running: http://127.0.0.1:3000");
  console.log("WORKDIR:", WORKDIR);
  console.log("QWEN_EXE:", QWEN_EXE);
  console.log("====================================\n");
});
