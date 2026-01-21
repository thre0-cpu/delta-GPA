import React, { useEffect, useRef } from 'react';
import ReactMarkdown from 'react-markdown';
import { Message } from '../types';
import { Bot, User, AlertCircle, Terminal } from 'lucide-react';

interface MessageListProps {
  messages: Message[];
  isGenerating: boolean;
}

const MessageList: React.FC<MessageListProps> = ({ messages, isGenerating }) => {
  const bottomRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages, messages.length, isGenerating]);

  if (messages.length === 0) {
    return (
      <div className="flex flex-col items-center justify-center h-full text-zinc-500 opacity-50">
        <Bot size={48} className="mb-4" />
        <p className="text-lg font-medium">Delta GPA</p>
        <p className="text-sm">Ready to assist.</p>
      </div>
    );
  }

  return (
    <div className="flex flex-col gap-6 p-4 md:p-6 max-w-3xl mx-auto w-full">
      {messages.map((msg) => {
        const isUser = msg.role === 'user';
        const isSystem = msg.role === 'system';

        if (isSystem) {
            return (
                <div key={msg.id} className="flex gap-3 justify-center text-xs text-zinc-500 my-2">
                    <span className="flex items-center gap-1 bg-zinc-900 border border-zinc-800 px-3 py-1 rounded-full">
                        {msg.content.includes("Error") ? <AlertCircle size={12} className="text-red-500"/> : <Terminal size={12} />}
                        {msg.content}
                    </span>
                </div>
            )
        }

        return (
          <div
            key={msg.id}
            className={`flex gap-4 ${isUser ? 'flex-row-reverse' : 'flex-row'}`}
          >
            <div
              className={`flex-shrink-0 w-8 h-8 rounded-full flex items-center justify-center ${
                isUser ? 'bg-zinc-700' : 'bg-primary'
              }`}
            >
              {isUser ? <User size={18} className="text-zinc-300" /> : <Bot size={18} className="text-white" />}
            </div>

            <div
              className={`flex flex-col max-w-[85%] md:max-w-[75%] ${
                isUser ? 'items-end' : 'items-start'
              }`}
            >
              <div
                className={`px-4 py-3 rounded-2xl shadow-sm overflow-hidden text-sm leading-relaxed ${
                  isUser
                    ? 'bg-surfaceHighlight text-zinc-100 rounded-tr-sm'
                    : 'bg-transparent text-zinc-100 px-0 py-0'
                }`}
              >
                {isUser ? (
                  <div className="whitespace-pre-wrap">{msg.content}</div>
                ) : (
                  <div className="prose prose-invert prose-sm max-w-none">
                     <ReactMarkdown
                        components={{
                            code({node, className, children, ...props}) {
                                const match = /language-(\w+)/.exec(className || '')
                                return match ? (
                                <div className="bg-zinc-950 rounded-md my-2 border border-zinc-800">
                                    <div className="flex justify-between px-3 py-1 text-xs text-zinc-500 bg-zinc-900 border-b border-zinc-800 rounded-t-md">
                                        <span>{match[1]}</span>
                                    </div>
                                    <pre className="p-3 overflow-x-auto m-0 !bg-transparent">
                                        <code className={className} {...props}>
                                            {children}
                                        </code>
                                    </pre>
                                </div>
                                ) : (
                                <code className="bg-zinc-800 text-zinc-200 px-1 py-0.5 rounded text-xs" {...props}>
                                    {children}
                                </code>
                                )
                            }
                        }}
                     >
                         {msg.content}
                     </ReactMarkdown>
                  </div>
                )}
              </div>
            </div>
          </div>
        );
      })}
      
      {isGenerating && messages[messages.length-1]?.role !== 'assistant' && (
          <div className="flex gap-4">
              <div className="flex-shrink-0 w-8 h-8 rounded-full bg-primary flex items-center justify-center">
                  <Bot size={18} className="text-white" />
              </div>
              <div className="flex items-center h-8">
                  <span className="w-2 h-2 bg-zinc-500 rounded-full animate-bounce mr-1"></span>
                  <span className="w-2 h-2 bg-zinc-500 rounded-full animate-bounce mr-1 delay-75"></span>
                  <span className="w-2 h-2 bg-zinc-500 rounded-full animate-bounce delay-150"></span>
              </div>
          </div>
      )}
      
      <div ref={bottomRef} />
    </div>
  );
};

export default MessageList;