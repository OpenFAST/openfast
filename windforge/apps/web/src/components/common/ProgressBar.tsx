import clsx from 'clsx';

type BarSize = 'sm' | 'md' | 'lg';

interface ProgressBarProps {
  value?: number;
  size?: BarSize;
  showLabel?: boolean;
  color?: 'accent' | 'green' | 'amber' | 'red' | 'blue';
  animated?: boolean;
  className?: string;
}

const SIZE_HEIGHTS: Record<BarSize, string> = {
  sm: 'h-1.5',
  md: 'h-3',
  lg: 'h-5',
};

const COLOR_BG: Record<string, string> = {
  accent: 'bg-accent-500',
  green: 'bg-emerald-500',
  amber: 'bg-amber-500',
  red: 'bg-red-500',
  blue: 'bg-blue-500',
};

const COLOR_GLOW: Record<string, string> = {
  accent: 'shadow-[0_0_8px_rgba(0,180,216,0.4)]',
  green: 'shadow-[0_0_8px_rgba(16,185,129,0.4)]',
  amber: 'shadow-[0_0_8px_rgba(245,158,11,0.4)]',
  red: 'shadow-[0_0_8px_rgba(239,68,68,0.4)]',
  blue: 'shadow-[0_0_8px_rgba(59,130,246,0.4)]',
};

export default function ProgressBar({
  value,
  size = 'md',
  showLabel = false,
  color = 'accent',
  animated = true,
  className,
}: ProgressBarProps) {
  const isIndeterminate = value === undefined || value === null;
  const clampedValue = isIndeterminate ? 0 : Math.min(100, Math.max(0, value));

  return (
    <div className={clsx('w-full', className)}>
      <div
        className={clsx(
          'w-full rounded-full bg-slate-700/50 overflow-hidden',
          SIZE_HEIGHTS[size],
        )}
      >
        {isIndeterminate ? (
          <div
            className={clsx(
              'h-full rounded-full',
              COLOR_BG[color],
              COLOR_GLOW[color],
              'animate-indeterminate',
            )}
            style={{ width: '40%' }}
          />
        ) : (
          <div
            className={clsx(
              'h-full rounded-full transition-all duration-500 ease-out',
              COLOR_BG[color],
              animated && clampedValue > 0 && clampedValue < 100 && COLOR_GLOW[color],
              animated && clampedValue > 0 && clampedValue < 100 && 'animate-pulse-slow',
            )}
            style={{ width: `${clampedValue}%` }}
          >
            {showLabel && size === 'lg' && clampedValue > 10 && (
              <span className="flex h-full items-center justify-end pr-2 text-[10px] font-semibold text-white">
                {Math.round(clampedValue)}%
              </span>
            )}
          </div>
        )}
      </div>
      {showLabel && size !== 'lg' && !isIndeterminate && (
        <div className="mt-1 text-right">
          <span className="text-xs font-medium text-slate-400">
            {Math.round(clampedValue)}%
          </span>
        </div>
      )}
    </div>
  );
}
