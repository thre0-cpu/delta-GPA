import React, { useState } from 'react';
import { ToolRequestPayload } from '../types';
import { X, ShieldAlert, Check, Play } from 'lucide-react';

interface ToolModalProps {
  request: ToolRequestPayload;
  onDecision: (approved: boolean, input: any) => void;
}

const ToolModal: React.FC<ToolModalProps> = ({ request, onDecision }) => {
  const [jsonInput, setJsonInput] = useState(JSON.stringify(request.input, null, 2));
  const [isEditing, setIsEditing] = useState(false);

  const handleApprove = () => {
    try {
      const parsed = JSON.parse(jsonInput);
      onDecision(true, parsed);
    } catch (e) {
      alert("Invalid JSON format in arguments.");
    }
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/70 backdrop-blur-sm p-4">
      <div className="bg-surface border border-zinc-700 rounded-xl shadow-2xl max-w-lg w-full overflow-hidden flex flex-col max-h-[90vh]">
        
        {/* Header */}
        <div className="flex items-center justify-between p-4 border-b border-zinc-700 bg-surfaceHighlight">
          <div className="flex items-center gap-2 text-amber-500">
            <ShieldAlert size={20} />
            <h3 className="font-semibold text-white">Tool Execution Request</h3>
          </div>
          <button onClick={() => onDecision(false, request.input)} className="text-zinc-400 hover:text-white transition-colors">
            <X size={20} />
          </button>
        </div>

        {/* Content */}
        <div className="p-6 overflow-y-auto">
          <p className="text-zinc-300 mb-4">
            The assistant wants to execute <span className="text-primary font-mono font-bold bg-blue-900/30 px-2 py-0.5 rounded">{request.toolName}</span>.
          </p>
          
          <div className="flex justify-between items-center mb-2">
            <span className="text-xs font-medium text-zinc-500 uppercase tracking-wider">Arguments</span>
            <button 
              onClick={() => setIsEditing(!isEditing)}
              className="text-xs text-primary hover:underline"
            >
              {isEditing ? "Reset" : "Edit JSON"}
            </button>
          </div>

          <div className="bg-zinc-950 rounded-lg border border-zinc-800 overflow-hidden">
            <textarea
              value={jsonInput}
              onChange={(e) => setJsonInput(e.target.value)}
              disabled={!isEditing}
              className={`w-full h-48 p-3 font-mono text-sm bg-transparent outline-none resize-none ${isEditing ? 'text-zinc-100' : 'text-zinc-400'}`}
              spellCheck={false}
            />
          </div>
          
          <p className="text-xs text-zinc-500 mt-2">
            Please review the arguments carefully before allowing execution.
          </p>
        </div>

        {/* Footer */}
        <div className="p-4 bg-zinc-900 border-t border-zinc-800 flex justify-end gap-3">
          <button
            onClick={() => onDecision(false, request.input)}
            className="px-4 py-2 rounded-lg text-sm font-medium text-zinc-300 hover:text-white hover:bg-zinc-800 transition-colors"
          >
            Deny
          </button>
          <button
            onClick={handleApprove}
            className="flex items-center gap-2 px-4 py-2 rounded-lg text-sm font-medium bg-primary text-white hover:bg-primaryHover transition-colors shadow-lg shadow-blue-900/20"
          >
            <Play size={16} fill="currentColor" />
            Allow Execution
          </button>
        </div>
      </div>
    </div>
  );
};

export default ToolModal;