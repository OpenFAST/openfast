import React, { useCallback, useRef } from 'react';
import Plot from 'react-plotly.js';
import { Download } from 'lucide-react';
import clsx from 'clsx';
import type { Data, Layout } from 'plotly.js';

interface PlotPanelProps {
  data: Data[];
  layout?: Partial<Layout>;
  title?: string;
  className?: string;
  height?: number | string;
  showExportButton?: boolean;
}

const BASE_LAYOUT: Partial<Layout> = {
  paper_bgcolor: 'transparent',
  plot_bgcolor: 'rgba(17,24,39,0.8)',
  font: {
    family: 'ui-monospace, monospace',
    size: 10,
    color: '#9ca3af',
  },
  margin: { t: 30, r: 20, b: 40, l: 50 },
  xaxis: {
    gridcolor: 'rgba(75,85,99,0.4)',
    zerolinecolor: 'rgba(75,85,99,0.6)',
    linecolor: 'rgba(75,85,99,0.6)',
  },
  yaxis: {
    gridcolor: 'rgba(75,85,99,0.4)',
    zerolinecolor: 'rgba(75,85,99,0.6)',
    linecolor: 'rgba(75,85,99,0.6)',
  },
  legend: {
    bgcolor: 'transparent',
    font: { size: 9, color: '#9ca3af' },
    orientation: 'h',
    y: -0.15,
    x: 0.5,
    xanchor: 'center',
  },
  showlegend: true,
};

export default function PlotPanel({
  data,
  layout,
  title,
  className,
  height = 250,
  showExportButton = true,
}: PlotPanelProps) {
  const plotRef = useRef<any>(null);

  const mergedLayout = React.useMemo<Partial<Layout>>(
    () => ({
      ...BASE_LAYOUT,
      ...layout,
      xaxis: { ...BASE_LAYOUT.xaxis, ...(layout?.xaxis || {}) },
      yaxis: { ...BASE_LAYOUT.yaxis, ...(layout?.yaxis || {}) },
      title: title
        ? {
            text: title,
            font: { size: 11, color: '#d1d5db' },
            x: 0.05,
            xanchor: 'left' as const,
          }
        : undefined,
    }),
    [layout, title]
  );

  const handleExport = useCallback(() => {
    const plotElement = plotRef.current?.el;
    if (!plotElement) return;

    import('plotly.js').then((Plotly) => {
      (Plotly as any).downloadImage(plotElement, {
        format: 'png',
        width: 800,
        height: 500,
        filename: title?.replace(/\s+/g, '_').toLowerCase() || 'plot',
      });
    });
  }, [title]);

  return (
    <div className={clsx('relative group', className)}>
      {showExportButton && (
        <button
          onClick={handleExport}
          className="absolute top-1 right-1 z-10 p-1.5 rounded
                     bg-gray-800/80 text-gray-400 hover:text-white
                     hover:bg-gray-700 opacity-0 group-hover:opacity-100
                     transition-all duration-200"
          title="Export as PNG"
        >
          <Download size={12} />
        </button>
      )}
      <Plot
        ref={plotRef}
        data={data}
        layout={mergedLayout}
        config={{
          displaylogo: false,
          responsive: true,
          modeBarButtonsToRemove: [
            'select2d',
            'lasso2d',
            'autoScale2d',
            'hoverClosestCartesian',
            'hoverCompareCartesian',
            'toggleSpikelines',
          ],
          displayModeBar: false,
        }}
        style={{ width: '100%', height: typeof height === 'number' ? `${height}px` : height }}
        useResizeHandler
      />
    </div>
  );
}
