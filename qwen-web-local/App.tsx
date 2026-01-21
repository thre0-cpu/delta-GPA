import React, { useState, KeyboardEvent, useEffect, useRef } from 'react';
import { Send, RefreshCw, Zap, ZapOff, Trash2, Power } from 'lucide-react';
import { useDeltaChat } from './hooks/useDeltaChat';
import MessageList from './components/MessageList';
import ToolModal from './components/ToolModal';

function App() {
  const { 
    sessionId, 
    messages, 
    status, 
    pendingTool, 
    toolEvents,
    sendMessage, 
    resetSession, 
    handleToolDecision,
    retryConnection
  } = useDeltaChat();

  const [input, setInput] = useState('');
  const [hasSentFirst, setHasSentFirst] = useState(false);
  const defaultPrompt = '你好，@sl/initialize.md 执行任务初始化';
  const defaultSummaryPrompt = '谢谢，@sl/summarize.md 总结任务';
  const toolEventsEndRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    toolEventsEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [toolEvents]);
  const handleSend = () => {
    if (!input.trim() || status === 'connecting') return;
    sendMessage(input);
    setInput('');
    setHasSentFirst(true);
  };

  const handleKeyDown = (e: KeyboardEvent<HTMLTextAreaElement>) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSend();
    }
  };

  const getStatusColor = () => {
    switch (status) {
      case 'connected': return 'text-green-500';
      case 'generating': return 'text-blue-500';
      case 'error': return 'text-red-500';
      case 'disconnected': return 'text-zinc-500';
      default: return 'text-yellow-500';
    }
  };

  return (
    <div className="flex flex-col h-screen bg-background text-zinc-100 font-sans overflow-hidden relative">
      
      {/* Top Bar */}
      <header className="flex-shrink-0 h-14 border-b border-zinc-800 flex items-center justify-between px-4 bg-surface/50 backdrop-blur-md z-10">
        <div className="flex items-center gap-3">
          <div className="w-8 h-8 bg-gradient-to-tr from-[#EE7788] to-[#f5a0b0] rounded-lg flex items-center justify-center shadow-lg shadow-blue-900/20">
            <span className="text-white text-sm tracking-tight" style={{ fontFamily: 'Cambria Math, STIX Two Text, serif' }}>𝓖</span>
          </div>
          <div>
            <h1 className="text-sm font-semibold tracking-wide">Delta GPA <span className="text-zinc-500 font-normal">V0.333</span></h1>
          </div>
        </div>

        <div className="flex items-center gap-3">
          <div className="hidden md:flex items-center gap-2 text-xs text-zinc-500 bg-zinc-900 px-3 py-1.5 rounded-full border border-zinc-800">
            <span className={`w-2 h-2 rounded-full ${status === 'connected' || status === 'generating' ? 'bg-green-500 animate-pulse' : 'bg-red-500'}`}></span>
            <span className="font-mono max-w-[80px] truncate" title={sessionId}>{sessionId}</span>
          </div>
          
          <button 
             onClick={retryConnection}
             title="Reconnect"
             className="p-2 text-zinc-400 hover:text-white hover:bg-zinc-800 rounded-md transition-all"
          >
             {status === 'disconnected' || status === 'error' ? <ZapOff size={18} /> : <Zap size={18} className={getStatusColor()}/>}
          </button>

          <button 
             onClick={resetSession}
             title="Reset Session"
             className="p-2 text-zinc-400 hover:text-red-400 hover:bg-zinc-800 rounded-md transition-all"
          >
            <Trash2 size={18} />
          </button>
        </div>
      </header>

      {/* Preset prompt banner */}
      <div className="px-4 pt-3">
        <div className="flex items-center justify-between border border-zinc-800 rounded-lg bg-surfaceHighlight/70 px-3 py-2 text-sm text-zinc-200">
          <span>任务填充</span>
          <div className="flex items-center gap-2">
            <button
              onClick={() => setInput(defaultPrompt)}
              disabled={hasSentFirst}
              className={`px-3 py-1 rounded-md text-white text-xs transition-all ${
                hasSentFirst ? 'opacity-50 cursor-not-allowed' : ''
              }`}
              style={{ 
                backgroundColor: hasSentFirst ? '#999999' : '#6677CC'
              }}
              onMouseEnter={(e) => {
                if (!hasSentFirst) e.currentTarget.style.backgroundColor = '#5566BB';
              }}
              onMouseLeave={(e) => {
                if (!hasSentFirst) e.currentTarget.style.backgroundColor = '#6677CC';
              }}
            >
              任务初始化
            </button>
            <button
              onClick={() => setInput(defaultSummaryPrompt)}
              disabled={!hasSentFirst}
              className={`px-3 py-1 rounded-md text-white text-xs transition-all ${
                !hasSentFirst ? 'opacity-50 cursor-not-allowed' : ''
              }`}
              style={{ 
                backgroundColor: !hasSentFirst ? '#999999' : '#EE7755'
              }}
              onMouseEnter={(e) => {
                if (hasSentFirst) e.currentTarget.style.backgroundColor = '#DD6644';
              }}
              onMouseLeave={(e) => {
                if (hasSentFirst) e.currentTarget.style.backgroundColor = '#EE7755';
              }}
            >
              任务总结
            </button>
          </div>
        </div>
      </div>

      {/* Main Chat Area */}
      <main className="flex-1 overflow-hidden grid grid-cols-1 lg:grid-cols-[2fr_1fr] gap-4 p-4">
        <div className="overflow-y-auto rounded-xl border border-zinc-800 bg-surface/40">
          <MessageList 
              messages={messages} 
              isGenerating={status === 'generating'} 
          />
        </div>

        <aside className="hidden lg:flex flex-col overflow-hidden rounded-xl border border-zinc-800 bg-surface/60">
          <div className="flex items-center justify-between px-4 py-3 border-b border-zinc-800">
            <div className="flex items-center gap-2 text-sm text-zinc-300">
              <Power size={16} className="text-primary" />
              <span>工具执行 / 任务</span>
            </div>
            <span className="text-xs text-zinc-500">{toolEvents.length} 条</span>
          </div>

          <div className="flex-1 overflow-y-auto px-4 py-3 space-y-3">
            {toolEvents.length === 0 && (
              <p className="text-xs text-zinc-500">暂无工具运行记录。</p>
            )}

            {toolEvents.slice(-50).map(evt => (
              <div key={evt.id} className="rounded-lg border border-zinc-800 bg-surfaceHighlight/60 p-3 text-xs text-zinc-200">
                <div className="flex items-center justify-between mb-1">
                  <span className={`font-semibold ${evt.kind === 'result' ? 'text-emerald-400' : 'text-blue-300'}`}>
                    {evt.kind === 'use' ? '工具调用' : (evt.isError ? '工具结果(错误)' : '工具结果')}
                  </span>
                  <span className="text-[10px] text-zinc-500 font-mono">{evt.toolUseId?.slice(0, 6) ?? 'n/a'}</span>
                </div>
                {evt.name && (
                  <div className="text-zinc-300 mb-2 font-mono text-[11px]">
                    <pre className="bg-zinc-950 p-2 rounded overflow-auto border border-zinc-800">{evt.name}</pre>
                  </div>
                )}
                {evt.input && (
                  <div className="text-zinc-400 whitespace-pre-wrap break-words mb-1 font-mono text-[10px]">
                    <span className="text-zinc-500">入参:</span>
                    {(() => {
                      try {
                        const parsed = typeof evt.input === 'string' ? JSON.parse(evt.input) : evt.input;
                        return <pre className="bg-zinc-950 p-2 rounded mt-1 overflow-auto">{JSON.stringify(parsed, null, 2)}</pre>;
                      } catch {
                        return <span>{typeof evt.input === 'string' ? evt.input : JSON.stringify(evt.input)}</span>;
                      }
                    })()}
                  </div>
                )}
                {evt.result && (
                  <div className="text-zinc-400 whitespace-pre-wrap break-words font-mono text-[10px]">
                    <span className="text-zinc-500">输出:</span>
                    {(() => {
                      try {
                        const parsed = typeof evt.result === 'string' ? JSON.parse(evt.result) : evt.result;
                        return <pre className="bg-zinc-950 p-2 rounded mt-1 overflow-auto">{JSON.stringify(parsed, null, 2)}</pre>;
                      } catch {
                        return <span>{typeof evt.result === 'string' ? evt.result : JSON.stringify(evt.result)}</span>;
                      }
                    })()}
                  </div>
                )}
              </div>
            ))}
            <div ref={toolEventsEndRef} />
          </div>
        </aside>
      </main>

      {/* Input Area */}
      <footer className="flex-shrink-0 p-4 border-t border-zinc-800 bg-surface/30 backdrop-blur-sm">
        <div className="max-w-3xl mx-auto relative">
          <div className="relative flex items-end gap-2 bg-surfaceHighlight rounded-xl border border-zinc-700 shadow-sm focus-within:ring-2 focus-within:ring-primary/50 focus-within:border-primary transition-all p-2">
            
            <textarea
              value={input}
              onChange={(e) => setInput(e.target.value)}
              onKeyDown={handleKeyDown}
              placeholder="Ask anything..."
              className="w-full bg-transparent text-sm text-zinc-100 placeholder-zinc-500 resize-none max-h-32 min-h-[44px] py-3 px-2 focus:outline-none scrollbar-hide"
              rows={1}
              style={{ height: 'auto', minHeight: '44px' }}
              // Auto-resize hack
              onInput={(e) => {
                const target = e.target as HTMLTextAreaElement;
                target.style.height = 'auto';
                target.style.height = `${Math.min(target.scrollHeight, 150)}px`;
              }}
            />
            
            <button
              onClick={handleSend}
              disabled={!input.trim() || status === 'generating' || status === 'connecting'}
              className={`p-2 rounded-lg mb-0.5 transition-all duration-200 ${
                input.trim() && status !== 'generating' 
                  ? 'bg-primary text-white shadow-md hover:bg-primaryHover' 
                  : 'bg-zinc-700 text-zinc-400 cursor-not-allowed'
              }`}
            >
              {status === 'generating' ? <RefreshCw size={18} className="animate-spin"/> : <Send size={18} />}
            </button>
          </div>
          
          <div className="text-center mt-2">
             <p className="text-[10px] text-zinc-600">
               Backend: <span className="font-mono">http://127.0.0.1:3000</span> • AI can make mistakes.
             </p>
          </div>
        </div>
      </footer>

      {/* Tool Approval Modal */}
      {pendingTool && (
        <ToolModal request={pendingTool} onDecision={handleToolDecision} />
      )}
    </div>
  );
}

export default App;