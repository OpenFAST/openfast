import clsx from 'clsx';

type BadgeSize = 'sm' | 'md' | 'lg';

interface StatusBadgeProps {
  status: string;
  size?: BadgeSize;
  className?: string;
}

const STATUS_COLORS: Record<string, { bg: string; text: string; ring: string }> = {
  pending: { bg: 'bg-slate-500/20', text: 'text-slate-300', ring: 'ring-slate-500/30' },
  queued: { bg: 'bg-blue-500/20', text: 'text-blue-300', ring: 'ring-blue-500/30' },
  running: { bg: 'bg-amber-500/20', text: 'text-amber-300', ring: 'ring-amber-500/30' },
  completed: { bg: 'bg-emerald-500/20', text: 'text-emerald-300', ring: 'ring-emerald-500/30' },
  failed: { bg: 'bg-red-500/20', text: 'text-red-300', ring: 'ring-red-500/30' },
  cancelled: { bg: 'bg-slate-500/20', text: 'text-slate-400', ring: 'ring-slate-500/30' },
  draft: { bg: 'bg-violet-500/20', text: 'text-violet-300', ring: 'ring-violet-500/30' },
  ready: { bg: 'bg-cyan-500/20', text: 'text-cyan-300', ring: 'ring-cyan-500/30' },
};

const SIZE_CLASSES: Record<BadgeSize, string> = {
  sm: 'px-1.5 py-0.5 text-[10px]',
  md: 'px-2.5 py-1 text-xs',
  lg: 'px-3 py-1.5 text-sm',
};

const DOT_SIZES: Record<BadgeSize, string> = {
  sm: 'h-1.5 w-1.5',
  md: 'h-2 w-2',
  lg: 'h-2.5 w-2.5',
};

export default function StatusBadge({ status, size = 'md', className }: StatusBadgeProps) {
  const normalizedStatus = status.toLowerCase();
  const colors = STATUS_COLORS[normalizedStatus] ?? STATUS_COLORS.pending;
  const isRunning = normalizedStatus === 'running';

  return (
    <span
      className={clsx(
        'inline-flex items-center gap-1.5 rounded-full font-medium ring-1 ring-inset capitalize',
        colors.bg,
        colors.text,
        colors.ring,
        SIZE_CLASSES[size],
        className,
      )}
    >
      <span
        className={clsx(
          'rounded-full',
          DOT_SIZES[size],
          isRunning ? 'animate-pulse bg-amber-400' : '',
          !isRunning && normalizedStatus === 'completed' && 'bg-emerald-400',
          !isRunning && normalizedStatus === 'failed' && 'bg-red-400',
          !isRunning && normalizedStatus === 'pending' && 'bg-slate-400',
          !isRunning && normalizedStatus === 'queued' && 'bg-blue-400',
          !isRunning && normalizedStatus === 'cancelled' && 'bg-slate-400',
          !isRunning && normalizedStatus === 'draft' && 'bg-violet-400',
          !isRunning && normalizedStatus === 'ready' && 'bg-cyan-400',
        )}
      />
      {status}
    </span>
  );
}
