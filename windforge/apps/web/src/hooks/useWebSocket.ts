import { useEffect, useRef, useState, useCallback } from 'react';
import type { WSMessage } from '@/types';

interface UseWebSocketOptions {
  simulationId: string | null;
  onMessage?: (message: WSMessage) => void;
  autoConnect?: boolean;
}

interface UseWebSocketReturn {
  isConnected: boolean;
  messages: WSMessage[];
  connect: () => void;
  disconnect: () => void;
  cancelCase: (caseId: string) => void;
}

export function useWebSocket({
  simulationId,
  onMessage,
  autoConnect = true,
}: UseWebSocketOptions): UseWebSocketReturn {
  const [isConnected, setIsConnected] = useState(false);
  const [messages, setMessages] = useState<WSMessage[]>([]);
  const wsRef = useRef<WebSocket | null>(null);
  const reconnectTimeoutRef = useRef<ReturnType<typeof setTimeout>>();
  const onMessageRef = useRef(onMessage);

  // Keep the callback ref up to date
  onMessageRef.current = onMessage;

  const connect = useCallback(() => {
    if (!simulationId) return;
    if (wsRef.current?.readyState === WebSocket.OPEN) return;

    // Determine WebSocket URL based on current location
    const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
    const host = window.location.host;
    const wsUrl = `${protocol}//${host}/ws/${simulationId}`;

    const ws = new WebSocket(wsUrl);
    wsRef.current = ws;

    ws.onopen = () => {
      setIsConnected(true);
      // Clear any pending reconnect
      if (reconnectTimeoutRef.current) {
        clearTimeout(reconnectTimeoutRef.current);
      }
    };

    ws.onmessage = (event) => {
      try {
        const message: WSMessage = JSON.parse(event.data);
        setMessages((prev) => [...prev, message]);
        onMessageRef.current?.(message);
      } catch {
        console.error('Failed to parse WebSocket message:', event.data);
      }
    };

    ws.onclose = (event) => {
      setIsConnected(false);
      wsRef.current = null;

      // Attempt to reconnect if not a clean close
      if (!event.wasClean && simulationId) {
        reconnectTimeoutRef.current = setTimeout(() => {
          connect();
        }, 3000);
      }
    };

    ws.onerror = () => {
      // onclose will be called after this
      console.error('WebSocket error');
    };
  }, [simulationId]);

  const disconnect = useCallback(() => {
    if (reconnectTimeoutRef.current) {
      clearTimeout(reconnectTimeoutRef.current);
    }
    if (wsRef.current) {
      wsRef.current.close(1000, 'Client disconnect');
      wsRef.current = null;
    }
    setIsConnected(false);
  }, []);

  const cancelCase = useCallback((caseId: string) => {
    if (wsRef.current?.readyState === WebSocket.OPEN) {
      wsRef.current.send(
        JSON.stringify({
          type: 'cancel_case',
          case_id: caseId,
        }),
      );
    }
  }, []);

  // Auto-connect when simulationId changes
  useEffect(() => {
    if (autoConnect && simulationId) {
      setMessages([]);
      connect();
    }

    return () => {
      disconnect();
    };
  }, [simulationId, autoConnect, connect, disconnect]);

  return {
    isConnected,
    messages,
    connect,
    disconnect,
    cancelCase,
  };
}
