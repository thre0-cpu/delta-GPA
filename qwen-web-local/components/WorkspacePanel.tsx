import React, { useState, useRef, useCallback } from 'react';
import { 
  FolderOpen, 
  RefreshCw, 
  Upload, 
  ChevronRight, 
  ChevronDown,
  File,
  FileText,
  FileCode,
  FileSpreadsheet,
  FileImage,
  Notebook,
  FolderClosed,
  AlertCircle
} from 'lucide-react';

interface FileNode {
  name: string;
  path: string;
  relativePath: string;
  isDirectory: boolean;
  children?: FileNode[];
}

interface WorkspacePanelProps {
  sessionId: string;
  workspacePath: string | null;
  files: FileNode[];
  isLoading: boolean;
  onRefresh: () => void;
  onOpenFile: (relativePath: string) => void;
  onUploadFiles: (files: FileList) => void;
}

// 根据文件扩展名返回对应图标
function getFileIcon(fileName: string) {
  const ext = fileName.split('.').pop()?.toLowerCase() || '';
  const iconClass = "flex-shrink-0";
  
  switch (ext) {
    case 'py':
      return <FileCode size={14} className={`${iconClass} text-yellow-400`} />;
    case 'ipynb':
      return <Notebook size={14} className={`${iconClass} text-orange-400`} />;
    case 'md':
    case 'txt':
      return <FileText size={14} className={`${iconClass} text-blue-300`} />;
    case 'csv':
    case 'tsv':
    case 'xlsx':
    case 'xls':
      return <FileSpreadsheet size={14} className={`${iconClass} text-green-400`} />;
    case 'json':
    case 'yaml':
    case 'yml':
    case 'toml':
      return <FileCode size={14} className={`${iconClass} text-purple-400`} />;
    case 'png':
    case 'jpg':
    case 'jpeg':
    case 'gif':
    case 'svg':
    case 'webp':
      return <FileImage size={14} className={`${iconClass} text-pink-400`} />;
    case 'html':
    case 'css':
    case 'js':
    case 'ts':
    case 'tsx':
    case 'jsx':
      return <FileCode size={14} className={`${iconClass} text-cyan-400`} />;
    default:
      return <File size={14} className={`${iconClass} text-zinc-400`} />;
  }
}

// 文件树节点组件
function FileTreeNode({ 
  node, 
  depth = 0, 
  onOpenFile 
}: { 
  node: FileNode; 
  depth?: number; 
  onOpenFile: (relativePath: string) => void;
}) {
  const [isExpanded, setIsExpanded] = useState(depth < 1); // 默认展开第一层
  
  const handleClick = () => {
    if (node.isDirectory) {
      setIsExpanded(!isExpanded);
    } else {
      onOpenFile(node.relativePath);
    }
  };
  
  return (
    <div>
      <div 
        className={`flex items-center gap-1.5 py-1 px-2 rounded cursor-pointer transition-colors hover:bg-zinc-700/50 ${
          !node.isDirectory ? 'hover:text-white' : ''
        }`}
        style={{ paddingLeft: `${depth * 12 + 8}px` }}
        onClick={handleClick}
        title={node.relativePath}
      >
        {node.isDirectory ? (
          <>
            {isExpanded ? (
              <ChevronDown size={12} className="flex-shrink-0 text-zinc-500" />
            ) : (
              <ChevronRight size={12} className="flex-shrink-0 text-zinc-500" />
            )}
            {isExpanded ? (
              <FolderOpen size={14} className="flex-shrink-0 text-amber-400" />
            ) : (
              <FolderClosed size={14} className="flex-shrink-0 text-amber-400" />
            )}
          </>
        ) : (
          <>
            <span className="w-3" /> {/* 占位 */}
            {getFileIcon(node.name)}
          </>
        )}
        <span className={`truncate text-xs ${node.isDirectory ? 'text-zinc-200' : 'text-zinc-300'}`}>
          {node.name}
        </span>
      </div>
      
      {node.isDirectory && isExpanded && node.children && (
        <div>
          {node.children.map((child, idx) => (
            <FileTreeNode 
              key={`${child.relativePath}-${idx}`} 
              node={child} 
              depth={depth + 1}
              onOpenFile={onOpenFile}
            />
          ))}
        </div>
      )}
    </div>
  );
}

export default function WorkspacePanel({
  sessionId,
  workspacePath,
  files,
  isLoading,
  onRefresh,
  onOpenFile,
  onUploadFiles
}: WorkspacePanelProps) {
  const fileInputRef = useRef<HTMLInputElement>(null);
  const [isDragging, setIsDragging] = useState(false);
  
  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    if (workspacePath) {
      setIsDragging(true);
    }
  }, [workspacePath]);
  
  const handleDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
  }, []);
  
  const handleDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDragging(false);
    
    if (workspacePath && e.dataTransfer.files.length > 0) {
      onUploadFiles(e.dataTransfer.files);
    }
  }, [workspacePath, onUploadFiles]);
  
  const handleFileSelect = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      onUploadFiles(e.target.files);
      e.target.value = ''; // 清空以便重复选择同一文件
    }
  };
  
  // 获取工作区文件夹名称
  const workspaceName = workspacePath?.split(/[/\\]/).pop() || '未设置';
  
  return (
    <div 
      className={`flex-1 rounded-xl border bg-surface/60 overflow-hidden flex flex-col min-h-0 transition-colors ${
        isDragging ? 'border-primary bg-primary/10' : 'border-zinc-800'
      }`}
      onDragOver={handleDragOver}
      onDragLeave={handleDragLeave}
      onDrop={handleDrop}
    >
      {/* 头部 */}
      <div className="flex items-center justify-between px-4 py-3 border-b border-zinc-800">
        <div className="flex items-center gap-2 text-sm text-zinc-300 min-w-0">
          <FolderOpen size={16} className="text-amber-400 flex-shrink-0" />
          <span className="truncate" title={workspacePath || '未设置工作区'}>
            {workspacePath ? `工作区（${workspaceName}）` : '工作区'}
          </span>
        </div>
        <div className="flex items-center gap-1">
          <button
            onClick={() => fileInputRef.current?.click()}
            disabled={!workspacePath}
            className={`p-1.5 rounded transition-colors ${
              workspacePath 
                ? 'text-zinc-400 hover:text-white hover:bg-zinc-700' 
                : 'text-zinc-600 cursor-not-allowed'
            }`}
            title="上传文件"
          >
            <Upload size={14} />
          </button>
          <button
            onClick={onRefresh}
            disabled={!workspacePath || isLoading}
            className={`p-1.5 rounded transition-colors ${
              workspacePath && !isLoading
                ? 'text-zinc-400 hover:text-white hover:bg-zinc-700' 
                : 'text-zinc-600 cursor-not-allowed'
            }`}
            title="刷新"
          >
            <RefreshCw size={14} className={isLoading ? 'animate-spin' : ''} />
          </button>
        </div>
        <input
          ref={fileInputRef}
          type="file"
          multiple
          className="hidden"
          onChange={handleFileSelect}
        />
      </div>
      
      {/* 文件列表 */}
      <div className="flex-1 overflow-y-auto py-2">
        {!workspacePath ? (
          <div className="flex flex-col items-center justify-center h-full text-zinc-500 px-4">
            <AlertCircle size={24} className="mb-2 text-zinc-600" />
            <p className="text-xs text-center">
              请先执行"任务初始化"<br/>创建工作区
            </p>
          </div>
        ) : isLoading ? (
          <div className="flex items-center justify-center h-full">
            <RefreshCw size={20} className="animate-spin text-zinc-500" />
          </div>
        ) : files.length === 0 ? (
          <div className="flex flex-col items-center justify-center h-full text-zinc-500 px-4">
            <FolderOpen size={24} className="mb-2 text-zinc-600" />
            <p className="text-xs text-center">
              工作区为空<br/>
              <span className="text-zinc-600">拖拽或点击上传文件</span>
            </p>
          </div>
        ) : (
          <div className="space-y-0.5">
            {files.map((node, idx) => (
              <FileTreeNode 
                key={`${node.relativePath}-${idx}`} 
                node={node} 
                onOpenFile={onOpenFile}
              />
            ))}
          </div>
        )}
        
        {/* 拖拽提示 */}
        {isDragging && (
          <div className="absolute inset-0 flex items-center justify-center bg-primary/20 border-2 border-dashed border-primary rounded-xl">
            <div className="text-center text-primary">
              <Upload size={32} className="mx-auto mb-2" />
              <p className="text-sm font-medium">释放以上传文件</p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
