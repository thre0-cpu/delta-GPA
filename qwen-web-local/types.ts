export type Role = 'user' | 'assistant' | 'system';

export interface Message {
  id: string;
  role: Role;
  content: string;
  timestamp: number;
}

export interface ToolRequestPayload {
  approvalId: string;
  toolName: string;
  input: any;
}

export type ToolEventKind = 'use' | 'result';

export interface ToolEvent {
  id: string;
  toolUseId?: string;
  name?: string;
  input?: any;
  result?: any;
  isError?: boolean;
  partial?: boolean;
  kind: ToolEventKind;
  timestamp: number;
}

export type ConnectionStatus = 'connecting' | 'connected' | 'generating' | 'done' | 'error' | 'disconnected';

export interface TodoItem {
  id: string;
  content: string;
  status: 'pending' | 'in-progress' | 'completed';
}

export interface ChatState {
  sessionId: string;
  status: ConnectionStatus;
  errorMessage: string | null;
}