import React, { useCallback, useMemo, useState } from 'react';
import { Copy, Check, FileText, Loader2 } from 'lucide-react';
import clsx from 'clsx';

interface FilePreviewProps {
  content: string | null;
  title?: string;
  isLoading?: boolean;
  className?: string;
}

function highlightLine(line: string): React.ReactNode {
  // Comment lines (start with ! or ---)
  if (/^\s*(!|-{3,})/.test(line)) {
    return <span className="text-green-400">{line}</span>;
  }
  // Section headers (lines with === or all caps section names)
  if (/^-{2,}\s+\S/.test(line) || /^[A-Z][A-Z_ ]{4,}/.test(line.trim())) {
    return <span className="text-blue-400 font-semibold">{line}</span>;
  }
  // Lines with key-value pairs (value then key)
  const kvMatch = line.match(/^(\s*)([\d.Ee+-]+)(\s+)(.+)$/);
  if (kvMatch) {
    return (
      <>
        <span>{kvMatch[1]}</span>
        <span className="text-amber-300">{kvMatch[2]}</span>
        <span>{kvMatch[3]}</span>
        <span className="text-gray-400">{kvMatch[4]}</span>
      </>
    );
  }
  // Boolean values
  if (/\b(True|False)\b/i.test(line)) {
    return (
      <span>
        {line.split(/\b(True|False)\b/i).map((part, i) =>
          /^(True|False)$/i.test(part) ? (
            <span key={i} className="text-purple-400">
              {part}
            </span>
          ) : (
            <span key={i}>{part}</span>
          )
        )}
      </span>
    );
  }
  return <span className="text-gray-300">{line}</span>;
}

export default function FilePreview({
  content,
  title,
  isLoading = false,
  className,
}: FilePreviewProps) {
  const [copied, setCopied] = useState(false);

  const lines = useMemo(() => {
    if (!content) return [];
    return content.split('\n');
  }, [content]);

  const handleCopy = useCallback(async () => {
    if (!content) return;
    try {
      await navigator.clipboard.writeText(content);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch {
      // Fallback for older browsers
      const textarea = document.createElement('textarea');
      textarea.value = content;
      document.body.appendChild(textarea);
      textarea.select();
      document.execCommand('copy');
      document.body.removeChild(textarea);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    }
  }, [content]);

  if (isLoading) {
    return (
      <div
        className={clsx(
          'bg-gray-900 border border-gray-700 rounded-lg overflow-hidden',
          className
        )}
      >
        {title && (
          <div className="flex items-center gap-2 px-4 py-2 bg-gray-800 border-b border-gray-700">
            <FileText size={14} className="text-gray-400" />
            <span className="text-xs font-medium text-gray-300">{title}</span>
          </div>
        )}
        <div className="p-4 space-y-2">
          <Loader2 size={16} className="animate-spin text-gray-500 mx-auto" />
          <div className="space-y-1.5">
            {Array.from({ length: 15 }).map((_, i) => (
              <div key={i} className="flex gap-3">
                <div
                  className="h-3 bg-gray-800 rounded animate-pulse"
                  style={{ width: '28px' }}
                />
                <div
                  className="h-3 bg-gray-800 rounded animate-pulse"
                  style={{ width: `${40 + Math.random() * 50}%` }}
                />
              </div>
            ))}
          </div>
        </div>
      </div>
    );
  }

  if (!content) {
    return (
      <div
        className={clsx(
          'bg-gray-900 border border-gray-700 rounded-lg overflow-hidden flex items-center justify-center p-8',
          className
        )}
      >
        <div className="text-center text-gray-500">
          <FileText size={32} className="mx-auto mb-2 opacity-40" />
          <p className="text-sm">No file content available</p>
          <p className="text-xs mt-1">Save the component to generate a preview</p>
        </div>
      </div>
    );
  }

  return (
    <div
      className={clsx(
        'bg-gray-900 border border-gray-700 rounded-lg overflow-hidden flex flex-col',
        className
      )}
    >
      {/* Header */}
      <div className="flex items-center gap-2 px-3 py-2 bg-gray-800 border-b border-gray-700 shrink-0">
        <FileText size={14} className="text-gray-400" />
        {title && (
          <span className="text-xs font-medium text-gray-300">{title}</span>
        )}
        <span className="ml-auto text-xs text-gray-500">
          {lines.length} lines
        </span>
        <button
          onClick={handleCopy}
          className="flex items-center gap-1 px-2 py-1 text-xs rounded
                     bg-gray-700 text-gray-300 hover:bg-gray-600
                     transition-colors"
          title="Copy to clipboard"
        >
          {copied ? (
            <>
              <Check size={12} className="text-green-400" />
              <span className="text-green-400">Copied</span>
            </>
          ) : (
            <>
              <Copy size={12} />
              <span>Copy</span>
            </>
          )}
        </button>
      </div>

      {/* Content */}
      <div className="overflow-auto flex-1">
        <pre className="text-xs leading-relaxed">
          <code>
            {lines.map((line, i) => (
              <div
                key={i}
                className="flex hover:bg-gray-800/50 transition-colors"
              >
                <span
                  className="w-10 shrink-0 text-right pr-3 text-gray-600
                             select-none border-r border-gray-800 font-mono"
                >
                  {i + 1}
                </span>
                <span className="pl-3 pr-4 font-mono whitespace-pre">
                  {highlightLine(line)}
                </span>
              </div>
            ))}
          </code>
        </pre>
      </div>
    </div>
  );
}
