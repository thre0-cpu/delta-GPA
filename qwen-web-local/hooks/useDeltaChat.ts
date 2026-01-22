import { useState, useEffect, useRef, useCallback } from 'react';
import { Message, ToolRequestPayload, ConnectionStatus, ToolEvent, TodoItem, FileNode } from '../types';

const API_BASE = "http://127.0.0.1:915";
const MAX_RECONNECT_ATTEMPTS = 3;
const RECONNECT_DELAYS = [500, 800, 1000]; // 快速重连：0.5s, 0.8s, 1s（共2.3秒）

export const useDeltaChat = () => {
  const [messages, setMessages] = useState<Message[]>([]);
  const [status, setStatus] = useState<ConnectionStatus>('connecting');
  const [sessionId, setSessionId] = useState<string>('');
  const [pendingTool, setPendingTool] = useState<ToolRequestPayload | null>(null);
  const [errorMsg, setErrorMsg] = useState<string | null>(null);
  const [toolEvents, setToolEvents] = useState<ToolEvent[]>([]);
  const [todoItems, setTodoItems] = useState<TodoItem[]>([]);
  const [connectionFailed, setConnectionFailed] = useState(false); // 连接彻底失败标志
  
  // 工作区相关状态
  const [workspacePath, setWorkspacePath] = useState<string | null>(null);
  const [workspaceFiles, setWorkspaceFiles] = useState<FileNode[]>([]);
  const [workspaceLoading, setWorkspaceLoading] = useState(false);

  const eventSourceRef = useRef<EventSource | null>(null);
  const reconnectAttemptsRef = useRef(0);
  const reconnectTimerRef = useRef<NodeJS.Timeout | null>(null);
  const connectSSERef = useRef<() => void>(() => {}); // 用于解决循环依赖
  const toolCallCounterRef = useRef(0); // 工具调用计数器
  const seenToolUseIdsRef = useRef<Set<string>>(new Set()); // 已处理过的 toolUseId
  const connectTimeoutRef = useRef<NodeJS.Timeout | null>(null); // 连接超时检测
  const refreshWorkspaceFilesRef = useRef<(() => void) | null>(null); // 用于在SSE回调中刷新文件列表

  // Initialize Session ID once
  useEffect(() => {
    const sid = crypto.randomUUID();
    setSessionId(sid);
  }, []);

  // 自动重连逻辑
  const attemptReconnect = useCallback(() => {
    if (reconnectAttemptsRef.current >= MAX_RECONNECT_ATTEMPTS) {
      console.log("[SSE] 重连次数已达上限，标记连接失败");
      setConnectionFailed(true);
      setStatus('disconnected');
      setErrorMsg(`连接失败：已尝试 ${MAX_RECONNECT_ATTEMPTS} 次重连均未成功`);
      return;
    }

    const delay = RECONNECT_DELAYS[reconnectAttemptsRef.current] || 3000;
    reconnectAttemptsRef.current++;
    
    console.log(`[SSE] 第 ${reconnectAttemptsRef.current} 次重连尝试，${delay}ms 后执行...`);
    setStatus('connecting');
    setErrorMsg(`正在重连... (第 ${reconnectAttemptsRef.current}/${MAX_RECONNECT_ATTEMPTS} 次)`);

    // 清除之前的定时器
    if (reconnectTimerRef.current) {
      clearTimeout(reconnectTimerRef.current);
    }

    reconnectTimerRef.current = setTimeout(() => {
      connectSSERef.current();
    }, delay);
  }, []);

  const connectSSE = useCallback(() => {
    if (!sessionId) return;
    
    // 清除之前的连接和超时
    if (eventSourceRef.current) {
      eventSourceRef.current.close();
      eventSourceRef.current = null;
    }
    if (connectTimeoutRef.current) {
      clearTimeout(connectTimeoutRef.current);
      connectTimeoutRef.current = null;
    }

    setStatus('connecting');

    const es = new EventSource(`${API_BASE}/api/sse?sessionId=${encodeURIComponent(sessionId)}`);
    eventSourceRef.current = es;

    // 连接超时检测：1.5秒内没收到 meta 事件就认为连接失败
    connectTimeoutRef.current = setTimeout(() => {
      if (es.readyState !== EventSource.OPEN || status === 'connecting') {
        console.log("[SSE] 连接超时，触发重连");
        es.close();
        eventSourceRef.current = null;
        attemptReconnect();
      }
    }, 1500);

    // 使用 onerror 属性
    es.onerror = () => {
      console.log("[SSE] onerror 触发, readyState:", es.readyState);
      // 立即关闭，防止 EventSource 自动重连
      es.close();
      eventSourceRef.current = null;
      
      // 清除超时检测（避免重复触发）
      if (connectTimeoutRef.current) {
        clearTimeout(connectTimeoutRef.current);
        connectTimeoutRef.current = null;
      }
      
      attemptReconnect();
    };

    es.addEventListener("meta", () => {
      // 连接成功，清除超时检测
      if (connectTimeoutRef.current) {
        clearTimeout(connectTimeoutRef.current);
        connectTimeoutRef.current = null;
      }
      
      console.log("[SSE] 收到 meta 事件，连接正常");
      reconnectAttemptsRef.current = 0;
      setConnectionFailed(false);
      setErrorMsg(null);
      setStatus("connected");
    });

    es.addEventListener("start", () => {
      setStatus("generating");
      // Prepare a new empty assistant message if the last one isn't already an empty assistant message
      setMessages(prev => {
        const last = prev[prev.length - 1];
        if (last && last.role === 'assistant' && last.content === '...') {
            return prev; // Already waiting
        }
        return [...prev, {
            id: crypto.randomUUID(),
            role: 'assistant',
            content: '',
            timestamp: Date.now()
        }];
      });
    });

    es.addEventListener("delta", (e) => {
      try {
        const { text } = JSON.parse(e.data);
        setMessages(prev => {
          const newArr = [...prev];
          const lastIdx = newArr.length - 1;
          if (lastIdx >= 0 && newArr[lastIdx].role === 'assistant') {
            newArr[lastIdx] = {
              ...newArr[lastIdx],
              content: newArr[lastIdx].content + text
            };
          }
          return newArr;
        });
      } catch (err) {
        console.error("Error parsing delta", err);
      }
    });

    es.addEventListener("done", () => {
      setStatus("connected");
      // 生成完成后自动刷新工作区文件列表
      setTimeout(() => {
        refreshWorkspaceFilesRef.current?.();
      }, 300);
    });

    es.addEventListener("interrupted", () => {
      setStatus("connected");
      // 添加系统消息提示用户已中止
      setMessages(prev => [...prev, {
        id: crypto.randomUUID(),
        role: 'system',
        content: '生成已中止',
        timestamp: Date.now()
      }]);
    });

    // 监听工作区创建事件（后端检测到工具创建 works/... 文件夹时推送）
    es.addEventListener("workspace", (e: MessageEvent) => {
      try {
        const { path: wsPath } = JSON.parse(e.data);
        console.log("[SSE] 收到 workspace 事件:", wsPath);
        if (wsPath) {
          setWorkspacePath(wsPath);
          // 延迟一点刷新文件列表，确保文件夹已创建
          setTimeout(() => {
            refreshWorkspaceFilesRef.current?.();
          }, 500);
        }
      } catch (err) {
        console.error("Error parsing workspace event", err);
      }
    });

    es.addEventListener("tool_event", (e: MessageEvent) => {
      try {
        const payload = JSON.parse(e.data);
        const toolUseId = payload.toolUseId;
        
        // 如果是 tool_use 事件且是首次收到该 toolUseId，在消息中插入标记
        if (payload.kind === 'use' && payload.name && toolUseId) {
          const isFirstSeen = !seenToolUseIdsRef.current.has(toolUseId);
          
          if (isFirstSeen) {
            seenToolUseIdsRef.current.add(toolUseId);
            toolCallCounterRef.current++;
            const toolIndex = toolCallCounterRef.current;
            // 使用 toolUseId 作为标记的唯一标识，方便后续查找
            const marker = `[[TOOL:${toolIndex}:${toolUseId}]]`;
            
            setMessages(prev => {
              const newArr = [...prev];
              const lastIdx = newArr.length - 1;
              if (lastIdx >= 0 && newArr[lastIdx].role === 'assistant') {
                newArr[lastIdx] = {
                  ...newArr[lastIdx],
                  content: newArr[lastIdx].content + marker
                };
              }
              return newArr;
            });
          }
        }
        
        // 更新或添加 toolEvent（用 toolUseId 去重更新）
        setToolEvents(prev => {
          // 查找是否已存在该 toolUseId 的事件
          const existingIdx = prev.findIndex(te => te.toolUseId === toolUseId && te.kind === payload.kind);
          
          if (existingIdx >= 0) {
            // 更新已存在的事件
            const updated = [...prev];
            updated[existingIdx] = {
              ...updated[existingIdx],
              input: payload.input ?? updated[existingIdx].input,
              result: payload.result ?? updated[existingIdx].result,
              isError: payload.isError ?? updated[existingIdx].isError,
              partial: payload.partial ?? false,
              timestamp: Date.now(),
            };
            return updated;
          } else {
            // 添加新事件
            return [...prev, {
              id: toolUseId || crypto.randomUUID(),
              toolUseId: toolUseId,
              name: payload.name,
              input: payload.input,
              result: payload.result,
              isError: payload.isError ?? false,
              partial: payload.partial ?? false,
              kind: payload.kind,
              timestamp: Date.now(),
            }];
          }
        });

        // 如果是 todo 相关工具调用，提取 todo 列表
        if (payload.name === 'todo_write' && payload.input) {
          try {
            const inputData = typeof payload.input === 'string' ? JSON.parse(payload.input) : payload.input;
            if (inputData.todos && Array.isArray(inputData.todos)) {
              setTodoItems(inputData.todos.map((t: any) => ({
                id: t.id || crypto.randomUUID(),
                content: t.content || t.title || '',
                status: t.status || 'pending',
              })));
            }
          } catch {
            // ignore parse error
          }
        }
        
        // 工具调用结果返回后刷新工作区（可能创建/修改了文件）
        if (payload.kind === 'result' && !payload.partial) {
          setTimeout(() => {
            refreshWorkspaceFilesRef.current?.();
          }, 300);
        }
      } catch (err) {
        console.error("Failed to parse tool_event", err);
      }
    });

    // 注意：error 处理已移至 es.onerror，此处仅处理自定义错误事件
    es.addEventListener("error", (e: Event) => {
       // 处理自定义的错误事件（后端发送的带 data 的错误）
       const messageEvent = e as MessageEvent;
       if (messageEvent.data) {
         try {
           const payload = JSON.parse(messageEvent.data);
           if (payload.message) {
             setErrorMsg(payload.message);
             setMessages(prev => [...prev, {
                id: crypto.randomUUID(),
                role: 'system',
                content: `Error: ${payload.message}`,
                timestamp: Date.now()
             }]);
           }
         } catch {
             // Ignore parse error
         }
       }
       // 重连逻辑已由 onerror 处理
    });

    es.addEventListener("tool_request", (e: MessageEvent) => {
      try {
        const payload = JSON.parse(e.data);
        setPendingTool(payload);
      } catch (err) {
        console.error("Failed to parse tool request", err);
      }
    });

    return () => {
      es.close();
    };
  }, [sessionId, attemptReconnect]);

  // 更新 connectSSE ref 以解决循环依赖
  useEffect(() => {
    connectSSERef.current = connectSSE;
  }, [connectSSE]);

  // Connect when sessionId is ready
  useEffect(() => {
    if (sessionId) {
      const cleanup = connectSSE();
      return () => {
        cleanup?.(); 
        if(eventSourceRef.current) eventSourceRef.current.close();
        // 清除重连定时器
        if(reconnectTimerRef.current) {
          clearTimeout(reconnectTimerRef.current);
          reconnectTimerRef.current = null;
        }
        // 清除连接超时检测
        if(connectTimeoutRef.current) {
          clearTimeout(connectTimeoutRef.current);
          connectTimeoutRef.current = null;
        }
      };
    }
  }, [sessionId, connectSSE]);

  const sendMessage = async (prompt: string, displayContent?: string) => {
    if (!prompt.trim() || !sessionId) return;

    // 检查 SSE 连接状态，如果断开则先重连
    if (!eventSourceRef.current || eventSourceRef.current.readyState === EventSource.CLOSED) {
      console.log("[sendMessage] SSE 连接已断开，尝试重连...");
      connectSSE();
      // 等待一小段时间让连接建立
      await new Promise(resolve => setTimeout(resolve, 500));
    }

    // Add User Message immediately (显示用户看到的内容，而非实际发送的内容)
    setMessages(prev => [...prev, {
      id: crypto.randomUUID(),
      role: 'user',
      content: displayContent ?? prompt,
      timestamp: Date.now()
    }]);

    setStatus('generating');

    try {
      const resp = await fetch(`${API_BASE}/api/chat`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sessionId, prompt }),
      });
      
      if (!resp.ok) {
        throw new Error(`Sending failed: ${resp.statusText}`);
      }
    } catch (err: any) {
      setErrorMsg(err.message);
      setStatus('error');
      
      // 发送失败可能是连接问题，尝试重连
      console.log("[sendMessage] 发送失败，尝试重连...");
      connectSSE();
      
      setMessages(prev => [...prev, {
        id: crypto.randomUUID(),
        role: 'system',
        content: `发送失败: ${err.message}，正在尝试重连...`,
        timestamp: Date.now()
      }]);
    }
  };

  const handleToolDecision = async (approved: boolean, inputData: any) => {
    if (!pendingTool) return;

    const { approvalId } = pendingTool;
    setPendingTool(null); // Close modal

    // Add a system note about the decision
    setMessages(prev => [...prev, {
        id: crypto.randomUUID(),
        role: 'system',
        content: approved 
            ? `Allowed execution of tool: **${pendingTool.toolName}**` 
            : `Denied execution of tool: **${pendingTool.toolName}**`,
        timestamp: Date.now()
    }]);

    try {
      await fetch(`${API_BASE}/api/tool/decision`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          approvalId,
          behavior: approved ? "allow" : "deny",
          updatedInput: inputData, // Pass back (potentially modified) input
          message: approved ? "" : "User denied via UI",
        }),
      });
    } catch (err) {
      console.error("Tool decision failed", err);
    }
  };

  const resetSession = async () => {
    setMessages([]);
    setErrorMsg(null);
    setToolEvents([]);
    setTodoItems([]);
    toolCallCounterRef.current = 0; // 重置工具调用计数器
    seenToolUseIdsRef.current.clear(); // 重置已处理的工具ID集合
    try {
      await fetch(`${API_BASE}/api/session/reset`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sessionId }),
      });
      setStatus('connected');
    } catch (err) {
      console.error("Reset failed", err);
    }
  };

  const retryConnection = () => {
      // 手动重试时重置计数和失败标志
      reconnectAttemptsRef.current = 0;
      setConnectionFailed(false);
      setErrorMsg(null);
      connectSSE();
  };

  const interruptGeneration = async () => {
    if (status !== 'generating') return;
    
    try {
      await fetch(`${API_BASE}/api/chat/interrupt`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sessionId }),
      });
    } catch (err) {
      console.error("Interrupt failed", err);
    }
  };

  // ========== 工作区功能 ==========
  
  // 设置工作区路径
  const setWorkspace = async (path: string) => {
    try {
      const resp = await fetch(`${API_BASE}/api/workspace/set`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sessionId, workspacePath: path }),
      });
      const data = await resp.json();
      if (data.ok) {
        setWorkspacePath(data.path);
        // 设置成功后自动刷新文件列表
        await refreshWorkspaceFiles();
      }
    } catch (err) {
      console.error("Set workspace failed", err);
    }
  };

  // 刷新工作区文件列表
  const refreshWorkspaceFiles = async () => {
    if (!sessionId) return;
    
    setWorkspaceLoading(true);
    try {
      const resp = await fetch(`${API_BASE}/api/workspace/files?sessionId=${encodeURIComponent(sessionId)}`);
      const data = await resp.json();
      if (data.ok) {
        setWorkspaceFiles(data.files || []);
        if (data.workspacePath) {
          setWorkspacePath(data.workspacePath);
        }
      }
    } catch (err) {
      console.error("Refresh workspace files failed", err);
    } finally {
      setWorkspaceLoading(false);
    }
  };

  // 更新 ref 以便在 SSE 回调中使用
  useEffect(() => {
    refreshWorkspaceFilesRef.current = refreshWorkspaceFiles;
  }, [sessionId]);

  // 打开工作区中的文件
  const openWorkspaceFile = async (relativePath: string) => {
    try {
      await fetch(`${API_BASE}/api/workspace/open`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sessionId, relativePath }),
      });
    } catch (err) {
      console.error("Open file failed", err);
    }
  };

  // 上传文件到工作区
  const uploadToWorkspace = async (files: FileList) => {
    if (!workspacePath) return;
    
    const formData = new FormData();
    formData.append('sessionId', sessionId);
    for (let i = 0; i < files.length; i++) {
      formData.append('files', files[i]);
    }
    
    try {
      const resp = await fetch(`${API_BASE}/api/workspace/upload`, {
        method: "POST",
        body: formData,
      });
      const data = await resp.json();
      if (data.ok) {
        // 上传成功后刷新文件列表
        await refreshWorkspaceFiles();
      }
    } catch (err) {
      console.error("Upload failed", err);
    }
  };

  // 注意：工作区路径现在通过 SSE 的 "workspace" 事件接收，
  // 后端在检测到工具创建 works/... 文件夹时自动推送

  return {
    sessionId,
    messages,
    status,
    errorMsg,
    pendingTool,
    toolEvents,
    todoItems,
    connectionFailed,
    // 工作区相关
    workspacePath,
    workspaceFiles,
    workspaceLoading,
    setWorkspace,
    refreshWorkspaceFiles,
    openWorkspaceFile,
    uploadToWorkspace,
    // 方法
    sendMessage,
    resetSession,
    handleToolDecision,
    retryConnection,
    interruptGeneration
  };
};