import { useState, useEffect, useRef, useCallback } from 'react';
import { Message, ToolRequestPayload, ConnectionStatus, ToolEvent } from '../types';

const API_BASE = "http://127.0.0.1:3000";

export const useDeltaChat = () => {
  const [messages, setMessages] = useState<Message[]>([]);
  const [status, setStatus] = useState<ConnectionStatus>('connecting');
  const [sessionId, setSessionId] = useState<string>('');
  const [pendingTool, setPendingTool] = useState<ToolRequestPayload | null>(null);
  const [errorMsg, setErrorMsg] = useState<string | null>(null);
  const [toolEvents, setToolEvents] = useState<ToolEvent[]>([]);

  const eventSourceRef = useRef<EventSource | null>(null);

  // Initialize Session ID once
  useEffect(() => {
    const sid = crypto.randomUUID();
    setSessionId(sid);
  }, []);

  const connectSSE = useCallback(() => {
    if (!sessionId) return;
    
    if (eventSourceRef.current) {
      eventSourceRef.current.close();
    }

    setStatus('connecting');
    setErrorMsg(null);

    const es = new EventSource(`${API_BASE}/api/sse?sessionId=${encodeURIComponent(sessionId)}`);
    eventSourceRef.current = es;

    es.addEventListener("meta", () => {
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
      setStatus("done");
    });

    es.addEventListener("tool_event", (e: MessageEvent) => {
      try {
        const payload = JSON.parse(e.data);
        setToolEvents(prev => [...prev, {
          id: crypto.randomUUID(),
          toolUseId: payload.toolUseId,
          name: payload.name,
          input: payload.input,
          result: payload.result,
          isError: payload.isError ?? false,
          partial: payload.partial ?? false,
          kind: payload.kind,
          timestamp: Date.now(),
        }]);
      } catch (err) {
        console.error("Failed to parse tool_event", err);
      }
    });

    es.addEventListener("error", (e: Event) => {
       // Check if it's a custom JSON error from backend or a generic network error
       // The original code tried to parse e.data, but standard EventSource error events don't usually have data property accessible this way in all browsers unless it's a custom event type named 'error'. 
       // However, to match original logic, we try to parse it.
       const messageEvent = e as MessageEvent;
       if (messageEvent.data) {
         try {
           const payload = JSON.parse(messageEvent.data);
           if (payload.message) {
             setErrorMsg(payload.message);
             // Add error as a system message
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
       
       if (es.readyState === EventSource.CLOSED) {
         setStatus("disconnected");
       } else {
         setStatus("error");
       }
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
  }, [sessionId]);

  // Connect when sessionId is ready
  useEffect(() => {
    if (sessionId) {
      const cleanup = connectSSE();
      return () => {
        cleanup?.(); 
        if(eventSourceRef.current) eventSourceRef.current.close();
      };
    }
  }, [sessionId, connectSSE]);

  const sendMessage = async (prompt: string) => {
    if (!prompt.trim() || !sessionId) return;

    // Add User Message immediately
    setMessages(prev => [...prev, {
      id: crypto.randomUUID(),
      role: 'user',
      content: prompt,
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
      setMessages(prev => [...prev, {
        id: crypto.randomUUID(),
        role: 'system',
        content: `Failed to send message: ${err.message}`,
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
      connectSSE();
  };

  return {
    sessionId,
    messages,
    status,
    errorMsg,
    pendingTool,
    toolEvents,
    sendMessage,
    resetSession,
    handleToolDecision,
    retryConnection
  };
};