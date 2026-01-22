import React, { useEffect, useRef, useState } from 'react';
import ReactMarkdown from 'react-markdown';
import { Message, ToolEvent } from '../types';
import { Bot, User, AlertCircle, Terminal, Wrench, X, FileText, AlertTriangle, Copy, Check } from 'lucide-react';

const API_BASE = "http://127.0.0.1:915";

// 检测消息是否可能是"幻觉执行"（声称完成但没有工具调用）
const detectHallucination = (content: string): boolean => {
  // 声称完成/执行的关键词
  const claimKeywords = [
    '已完成', '已创建', '已生成', '已执行', '已运行', '已保存', '已写入',
    '创建了', '生成了', '执行了', '运行了', '保存了', '写入了',
    '完成了', '建立了', '添加了', '修改了', '更新了', '删除了',
    '我创建', '我生成', '我执行', '我运行', '我保存', '我写入',
    '文件已', '代码已', '脚本已', 'Notebook已',
  ];
  
  // 检查是否有工具调用标记
  const hasToolMarker = /\[\[TOOL:\d+:[^\]]+\]\]/.test(content);
  
  // 检查是否有声称完成的关键词
  const hasClaim = claimKeywords.some(kw => content.includes(kw));
  
  // 如果有声称但没有工具调用，可能是幻觉
  return hasClaim && !hasToolMarker;
};

// 检测字符串是否是文件路径（包括纯文件名）
const isFilePath = (text: string): boolean => {
  const trimmed = text.trim();
  
  // 匹配常见文件路径模式
  const patterns = [
    /^[A-Za-z]:[\\/]/,  // Windows 绝对路径: D:\xxx\yyy.ext
    /^\.{0,2}[\\/]/,     // 相对路径 ./ 或 ../
    /^works[\\/]/,       // works 目录下的文件
    /^prompts[\\/]/,     // prompts 目录下的文件
    /^APIs[\\/]/,        // APIs 目录下的文件
    /^envs[\\/]/,        // envs 目录下的文件
    /\.(py|ipynb|md|txt|csv|json|ts|tsx|js|jsx|html|css|yaml|yml|toml|sh|bat)$/i,  // 常见文件扩展名（包括纯文件名）
  ];
  return patterns.some(p => p.test(trimmed));
};

// 从消息列表中提取最新的工作文件夹路径
const extractCurrentWorkDir = (messages: Message[]): string => {
  // 从后往前遍历消息，找到最近提到的工作文件夹
  for (let i = messages.length - 1; i >= 0; i--) {
    const content = messages[i].content;
    // 匹配多种格式的工作目录路径：
    // - works/20260122_005_newtask
    // - works\20260122_005_newtask
    // - 20260122_005_newtask（纯目录名）
    const patterns = [
      /works[\\\/](\d{8}_\d{3}_[a-zA-Z]+)/g,  // 完整路径
      /(\d{8}_\d{3}_newtask)/g,                // 纯目录名
    ];
    
    for (const pattern of patterns) {
      const matches = content.match(pattern);
      if (matches && matches.length > 0) {
        let result = matches[matches.length - 1];
        // 确保返回完整的 works/xxx 格式
        if (!result.startsWith('works')) {
          result = 'works/' + result;
        }
        console.log("[extractCurrentWorkDir] 找到工作目录:", result);
        return result.replace(/\\/g, '/');
      }
    }
  }
  console.log("[extractCurrentWorkDir] 未找到工作目录");
  return '';
};

// 打开本地文件（支持上下文目录）
const openLocalFile = async (filePath: string, contextDir?: string) => {
  try {
    const response = await fetch(`${API_BASE}/api/open-file`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ filePath, contextDir }),
    });
    const result = await response.json();
    if (!response.ok) {
      console.error('打开文件失败:', result.error);
      alert(`无法打开文件: ${result.error}\n路径: ${result.path || filePath}`);
    }
  } catch (err) {
    console.error('打开文件请求失败:', err);
  }
};

interface MessageListProps {
  messages: Message[];
  isGenerating: boolean;
  toolEvents?: ToolEvent[];
}

// 工具引用弹窗组件
interface ToolRefPopupProps {
  toolEvent: ToolEvent | null;
  toolIndex: number;
  position: { x: number; y: number };
  onClose: () => void;
}

const ToolRefPopup: React.FC<ToolRefPopupProps> = ({ toolEvent, toolIndex, position, onClose }) => {
  const [copiedField, setCopiedField] = useState<'name' | 'input' | 'result' | null>(null);

  if (!toolEvent) return null;

  const handleCopy = async (text: string, field: 'name' | 'input' | 'result') => {
    try {
      await navigator.clipboard.writeText(text);
      setCopiedField(field);
      setTimeout(() => setCopiedField(null), 1500);
    } catch (err) {
      console.error('复制失败:', err);
    }
  };

  const formatInput = (input: any): string => {
    if (typeof input === 'string') {
      try {
        return JSON.stringify(JSON.parse(input), null, 2);
      } catch {
        return input;
      }
    }
    return JSON.stringify(input, null, 2);
  };

  const formatResult = (result: any): string => {
    if (typeof result === 'string') return result;
    return JSON.stringify(result, null, 2);
  };

  return (
    <div 
      className="fixed z-50 bg-zinc-900 border border-zinc-700 rounded-lg shadow-xl max-w-md"
      style={{ 
        left: Math.min(position.x, window.innerWidth - 420), 
        top: Math.min(position.y + 10, window.innerHeight - 300)
      }}
      onClick={(e) => e.stopPropagation()}
    >
      <div className="flex items-center justify-between px-3 py-2 border-b border-zinc-700 bg-zinc-800 rounded-t-lg">
        <div className="flex items-center gap-2 text-sm font-medium text-zinc-200">
          <Wrench size={14} className="text-blue-400" />
          <span>工具调用 #{toolIndex}</span>
        </div>
        <button onClick={onClose} className="text-zinc-400 hover:text-white">
          <X size={14} />
        </button>
      </div>
      <div className="p-3 space-y-2 max-h-64 overflow-auto">
        <div className="flex items-center justify-between group">
          <div>
            <span className="text-xs text-zinc-500">工具名称：</span>
            <span className="text-sm text-blue-400 font-mono ml-1">{toolEvent.name || '未知'}</span>
          </div>
          <button 
            onClick={() => handleCopy(toolEvent.name || '', 'name')}
            className="opacity-0 group-hover:opacity-100 p-1 text-zinc-400 hover:text-white hover:bg-zinc-700 rounded transition-all"
            title="复制工具名称"
          >
            {copiedField === 'name' ? <Check size={12} className="text-green-400" /> : <Copy size={12} />}
          </button>
        </div>
        {toolEvent.input && (
          <div>
            <div className="flex items-center justify-between group">
              <span className="text-xs text-zinc-500">入参：</span>
              <button 
                onClick={() => handleCopy(formatInput(toolEvent.input), 'input')}
                className="opacity-0 group-hover:opacity-100 p-1 text-zinc-400 hover:text-white hover:bg-zinc-700 rounded transition-all"
                title="复制入参"
              >
                {copiedField === 'input' ? <Check size={12} className="text-green-400" /> : <Copy size={12} />}
              </button>
            </div>
            <pre className="mt-1 text-xs bg-zinc-950 p-2 rounded overflow-auto max-h-32 text-zinc-300 cursor-pointer hover:bg-zinc-900 transition-colors"
              onClick={() => handleCopy(formatInput(toolEvent.input), 'input')}
              title="点击复制"
            >
              {formatInput(toolEvent.input)}
            </pre>
          </div>
        )}
        {toolEvent.result && (
          <div>
            <div className="flex items-center justify-between group">
              <span className="text-xs text-zinc-500">结果：</span>
              <button 
                onClick={() => handleCopy(formatResult(toolEvent.result), 'result')}
                className="opacity-0 group-hover:opacity-100 p-1 text-zinc-400 hover:text-white hover:bg-zinc-700 rounded transition-all"
                title="复制结果"
              >
                {copiedField === 'result' ? <Check size={12} className="text-green-400" /> : <Copy size={12} />}
              </button>
            </div>
            <pre className="mt-1 text-xs bg-zinc-950 p-2 rounded overflow-auto max-h-32 text-zinc-300 cursor-pointer hover:bg-zinc-900 transition-colors"
              onClick={() => handleCopy(formatResult(toolEvent.result), 'result')}
              title="点击复制"
            >
              {formatResult(toolEvent.result).slice(0, 500) + (formatResult(toolEvent.result).length > 500 ? '...' : '')}
            </pre>
          </div>
        )}
      </div>
    </div>
  );
};

// 解析消息内容，将工具引用标记转换为可点击元素
const parseToolReferences = (
  content: string, 
  toolEvents: ToolEvent[], 
  onRefClick: (toolEvent: ToolEvent, index: number, e: React.MouseEvent) => void
) => {
  const regex = /\[\[TOOL:(\d+):([^\]]+)\]\]/g;
  const parts: (string | React.ReactNode)[] = [];
  let lastIndex = 0;
  let match;

  while ((match = regex.exec(content)) !== null) {
    // 添加标记前的文本
    if (match.index > lastIndex) {
      parts.push(content.slice(lastIndex, match.index));
    }
    
    const toolIndex = parseInt(match[1], 10);
    const toolUseId = match[2];
    // 用 toolUseId 查找对应的工具事件（只找 kind='use' 的）
    const toolEvent = toolEvents.find(te => te.toolUseId === toolUseId && te.kind === 'use');
    
    // 添加可点击的工具引用标识
    parts.push(
      <sup
        key={toolUseId}
        className="inline-flex items-center justify-center w-4 h-4 ml-0.5 text-[10px] font-bold text-white bg-blue-500 hover:bg-blue-400 rounded cursor-pointer transition-colors"
        title={`工具调用 #${toolIndex}: ${toolEvent?.name || '未知'}`}
        onClick={(e) => {
          e.stopPropagation();
          if (toolEvent) onRefClick(toolEvent, toolIndex, e);
        }}
      >
        {toolIndex}
      </sup>
    );
    
    lastIndex = match.index + match[0].length;
  }
  
  // 添加剩余文本
  if (lastIndex < content.length) {
    parts.push(content.slice(lastIndex));
  }
  
  return parts.length > 0 ? parts : [content];
};

const MessageList: React.FC<MessageListProps> = ({ messages, isGenerating, toolEvents = [] }) => {
  const bottomRef = useRef<HTMLDivElement>(null);
  const [popupInfo, setPopupInfo] = useState<{
    toolEvent: ToolEvent;
    toolIndex: number;
    position: { x: number; y: number };
  } | null>(null);

  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages, messages.length, isGenerating]);

  const handleRefClick = (toolEvent: ToolEvent, toolIndex: number, e: React.MouseEvent) => {
    setPopupInfo({
      toolEvent,
      toolIndex,
      position: { x: e.clientX, y: e.clientY }
    });
  };

  const closePopup = () => setPopupInfo(null);

  // 点击其他地方关闭弹窗
  useEffect(() => {
    const handleClickOutside = () => {
      if (popupInfo) closePopup();
    };
    document.addEventListener('click', handleClickOutside);
    return () => document.removeEventListener('click', handleClickOutside);
  }, [popupInfo]);

  if (messages.length === 0) {
    return (
      <div className="flex flex-col items-center justify-center h-full text-zinc-500 opacity-50">
        <Bot size={48} className="mb-4" />
        <p className="text-lg font-medium">Delta GPA</p>
        <p className="text-sm">Ready to assist.</p>
      </div>
    );
  }

  // 获取当前工作文件夹
  const currentWorkDir = extractCurrentWorkDir(messages);

  // 自定义渲染器，支持工具引用
  const renderContentWithRefs = (content: string) => {
    // 先移除工具标记用于 Markdown 渲染，但保留位置信息
    const cleanContent = content.replace(/\[\[TOOL:\d+:[^\]]+\]\]/g, '');
    const hasToolRefs = content !== cleanContent;
    
    // 自定义 code 组件，支持文件路径点击
    const codeComponent = ({node, className, children, ...props}: any) => {
      const match = /language-(\w+)/.exec(className || '');
      const codeText = String(children).replace(/\n$/, '');
      
      // 代码块
      if (match) {
        return (
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
        );
      }
      
      // 内联代码 - 检查是否是文件路径
      if (isFilePath(codeText)) {
        return (
          <code 
            className="bg-zinc-800 text-blue-400 px-1 py-0.5 rounded text-xs cursor-pointer hover:bg-zinc-700 hover:text-blue-300 transition-colors inline-flex items-center gap-1"
            onClick={() => openLocalFile(codeText, currentWorkDir)}
            title={`点击在 VS Code 中打开${currentWorkDir ? ` (工作目录: ${currentWorkDir})` : ''}`}
            {...props}
          >
            <FileText size={12} className="inline" />
            {children}
          </code>
        );
      }
      
      // 普通内联代码
      return (
        <code className="bg-zinc-800 text-zinc-200 px-1 py-0.5 rounded text-xs" {...props}>
          {children}
        </code>
      );
    };
    
    if (!hasToolRefs) {
      return (
        <ReactMarkdown components={{ code: codeComponent }}>
          {content}
        </ReactMarkdown>
      );
    }

    // 有工具引用时，用自定义解析
    const parts = parseToolReferences(content, toolEvents, handleRefClick);
    
    return (
      <>
        {parts.map((part, idx) => {
          if (typeof part === 'string') {
            return (
              <ReactMarkdown key={idx} components={{ code: codeComponent }}>
                {part}
              </ReactMarkdown>
            );
          }
          return part;
        })}
      </>
    );
  };

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
                className={`px-4 py-3 rounded-2xl shadow-sm text-sm leading-relaxed ${
                  isUser
                    ? 'bg-surfaceHighlight text-zinc-100 rounded-tr-sm'
                    : 'bg-zinc-900/80 text-zinc-100 border border-zinc-700 rounded-tl-sm'
                }`}
              >
                {isUser ? (
                  <div className="whitespace-pre-wrap">{msg.content}</div>
                ) : (
                  <div className="prose prose-invert prose-sm max-w-none">
                    {renderContentWithRefs(msg.content)}
                  </div>
                )}
              </div>
              
              {/* 幻觉警告：声称完成但没有工具调用 */}
              {!isUser && detectHallucination(msg.content) && (
                <div className="flex items-center gap-2 mt-2 px-3 py-1.5 bg-amber-900/30 border border-amber-600/50 rounded-lg text-xs text-amber-400">
                  <AlertTriangle size={14} className="flex-shrink-0" />
                  <span>⚠️ 此消息声称执行了操作，但未检测到工具调用，可能是幻觉输出</span>
                </div>
              )}
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
      
      {/* 工具引用弹窗 */}
      {popupInfo && (
        <ToolRefPopup
          toolEvent={popupInfo.toolEvent}
          toolIndex={popupInfo.toolIndex}
          position={popupInfo.position}
          onClose={closePopup}
        />
      )}
    </div>
  );
};

export default MessageList;