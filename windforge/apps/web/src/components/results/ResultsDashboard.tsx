import { useState, useEffect, useMemo, useCallback } from 'react';
import { useParams } from 'react-router-dom';
import {
  BarChart3,
  AlertTriangle,
  TrendingUp,
  Layers,
  FileText,
  Download,
  ChevronDown,
  Maximize2,
  Check,
} from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';
import Plot from 'react-plotly.js';
import apiClient from '@/api/client';

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface ChannelStatistic {
  min: number;
  max: number;
  mean: number;
  std: number;
  abs_max?: number;
}

interface ResultsStatistics {
  id: string;
  simulation_case_id: string;
  dlc_number: string;
  wind_speed: number;
  channel_statistics: Record<string, ChannelStatistic> | null;
}

interface ResultsDEL {
  id: string;
  simulation_id: string;
  del_results: Record<string, number> | null;
  m_exponent: number;
  n_equivalent: number;
}

interface ResultsExtreme {
  id: string;
  simulation_id: string;
  extreme_loads: Record<
    string,
    { max: number; min: number; safety_factor: number; design_value: number }
  > | null;
}

interface Simulation {
  id: string;
  name: string;
  status: string;
}

// ---------------------------------------------------------------------------
// Tabs
// ---------------------------------------------------------------------------

type TabKey = 'statistics' | 'del' | 'extreme' | 'envelopes' | 'report';

const TABS: { key: TabKey; label: string; icon: React.ReactNode }[] = [
  { key: 'statistics', label: 'Statistics', icon: <BarChart3 size={14} /> },
  { key: 'del', label: 'DEL', icon: <TrendingUp size={14} /> },
  { key: 'extreme', label: 'Extreme Loads', icon: <AlertTriangle size={14} /> },
  { key: 'envelopes', label: 'Load Envelopes', icon: <Layers size={14} /> },
  { key: 'report', label: 'Report', icon: <FileText size={14} /> },
];

// ---------------------------------------------------------------------------
// Commonly referenced channels
// ---------------------------------------------------------------------------

const COMMON_CHANNELS = [
  'GenPwr',
  'RotSpeed',
  'BldPitch1',
  'TwrBsMxt',
  'TwrBsMyt',
  'RootMxc1',
  'RootMyc1',
  'YawBrMxp',
  'YawBrMyp',
  'LSShftFxa',
  'LSShftTq',
  'TTDspFA',
  'TTDspSS',
];

const DLC_COLORS: Record<string, string> = {
  '1.1': '#00b4d8',
  '1.2': '#0088a8',
  '1.3': '#ef4444',
  '1.4': '#f97316',
  '1.5': '#eab308',
  '2.1': '#8b5cf6',
  '2.2': '#a855f7',
  '2.3': '#d946ef',
  '3.1': '#10b981',
  '3.2': '#14b8a6',
  '3.3': '#06b6d4',
  '4.1': '#84cc16',
  '4.2': '#22c55e',
  '5.1': '#f43f5e',
  '6.1': '#6366f1',
  '6.2': '#818cf8',
  '6.3': '#a5b4fc',
  '6.4': '#c4b5fd',
  '7.1': '#fb923c',
  '8.1': '#fbbf24',
};

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------

export default function ResultsDashboard() {
  const { projectId } = useParams<{ projectId: string }>();

  const [activeTab, setActiveTab] = useState<TabKey>('statistics');
  const [simulations, setSimulations] = useState<Simulation[]>([]);
  const [selectedSimId, setSelectedSimId] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  // Data
  const [statistics, setStatistics] = useState<ResultsStatistics[]>([]);
  const [delResults, setDelResults] = useState<ResultsDEL | null>(null);
  const [extremeResults, setExtremeResults] = useState<ResultsExtreme | null>(null);

  // Filters
  const [selectedChannel, setSelectedChannel] = useState('TwrBsMxt');
  const [delMExponent, setDelMExponent] = useState(4);

  // Report config
  const [reportConfig, setReportConfig] = useState({
    includeDEL: true,
    includeExtreme: true,
    includeEnvelopes: true,
  });
  const [generatingReport, setGeneratingReport] = useState(false);

  // Envelope axes
  const [envelopeXChannel, setEnvelopeXChannel] = useState('TwrBsMxt');
  const [envelopeYChannel, setEnvelopeYChannel] = useState('TwrBsMyt');

  // ---- Fetch simulations ----
  useEffect(() => {
    if (!projectId) return;
    setLoading(true);
    apiClient
      .get(`/projects/${projectId}/simulations`)
      .then((res) => {
        const sims = (res.data as Simulation[]).filter((s) => s.status === 'completed');
        setSimulations(sims);
        if (sims.length > 0) setSelectedSimId(sims[0].id);
      })
      .catch(() => toast.error('Failed to load simulations'))
      .finally(() => setLoading(false));
  }, [projectId]);

  // ---- Fetch results when sim changes ----
  useEffect(() => {
    if (!projectId || !selectedSimId) return;

    apiClient
      .get(`/projects/${projectId}/simulations/${selectedSimId}/results/statistics`)
      .then((res) => setStatistics(res.data as ResultsStatistics[]))
      .catch(() => setStatistics([]));

    apiClient
      .get(`/projects/${projectId}/simulations/${selectedSimId}/results/del`)
      .then((res) => setDelResults(res.data as ResultsDEL))
      .catch(() => setDelResults(null));

    apiClient
      .get(`/projects/${projectId}/simulations/${selectedSimId}/results/extreme`)
      .then((res) => setExtremeResults(res.data as ResultsExtreme))
      .catch(() => setExtremeResults(null));
  }, [projectId, selectedSimId]);

  // ---- All channels from statistics ----
  const allChannels = useMemo(() => {
    const channelSet = new Set<string>();
    statistics.forEach((s) => {
      if (s.channel_statistics) {
        Object.keys(s.channel_statistics).forEach((ch) => channelSet.add(ch));
      }
    });
    return Array.from(channelSet).sort();
  }, [statistics]);

  // ---- Statistics for selected channel ----
  const channelStats = useMemo(() => {
    return statistics
      .filter((s) => s.channel_statistics?.[selectedChannel])
      .map((s) => ({
        dlc: s.dlc_number,
        windSpeed: s.wind_speed,
        ...s.channel_statistics![selectedChannel],
      }))
      .sort((a, b) => {
        if (a.dlc !== b.dlc) return a.dlc.localeCompare(b.dlc);
        return a.windSpeed - b.windSpeed;
      });
  }, [statistics, selectedChannel]);

  // ---- Aggregated statistics ----
  const aggregatedStats = useMemo(() => {
    if (channelStats.length === 0) return null;
    return {
      min: Math.min(...channelStats.map((s) => s.min)),
      max: Math.max(...channelStats.map((s) => s.max)),
      mean: channelStats.reduce((sum, s) => sum + s.mean, 0) / channelStats.length,
      std: Math.sqrt(
        channelStats.reduce((sum, s) => sum + s.std * s.std, 0) / channelStats.length,
      ),
      abs_max: Math.max(...channelStats.map((s) => s.abs_max ?? Math.max(Math.abs(s.min), Math.abs(s.max)))),
    };
  }, [channelStats]);

  // ---- Box plot data ----
  const boxPlotData = useMemo(() => {
    const windSpeedBins = [...new Set(channelStats.map((s) => s.windSpeed))].sort(
      (a, b) => a - b,
    );

    return windSpeedBins.map((ws) => {
      const vals = channelStats
        .filter((s) => s.windSpeed === ws)
        .map((s) => s.max);
      return {
        y: vals,
        type: 'box' as const,
        name: `${ws} m/s`,
        marker: { color: '#00b4d8' },
        boxpoints: 'all' as const,
        jitter: 0.3,
        pointpos: -1.5,
      };
    });
  }, [channelStats]);

  // ---- DEL bar chart ----
  const delBarData = useMemo(() => {
    if (!delResults?.del_results) return [];
    const entries = Object.entries(delResults.del_results).sort((a, b) => b[1] - a[1]);
    return [
      {
        x: entries.map(([ch]) => ch),
        y: entries.map(([, val]) => val),
        type: 'bar' as const,
        marker: {
          color: entries.map(([ch]) => DLC_COLORS[ch] ?? '#00b4d8'),
        },
      },
    ];
  }, [delResults]);

  // ---- Extreme loads sorted by utilization ----
  const extremeEntries = useMemo(() => {
    if (!extremeResults?.extreme_loads) return [];
    return Object.entries(extremeResults.extreme_loads)
      .map(([channel, data]) => ({
        channel,
        ...data,
        utilization: data.design_value > 0 ? data.max / data.design_value : 0,
      }))
      .sort((a, b) => b.utilization - a.utilization);
  }, [extremeResults]);

  // ---- Generate scatter plot for envelopes ----
  const envelopeTraces = useMemo(() => {
    const dlcGroups: Record<string, { x: number[]; y: number[] }> = {};

    statistics.forEach((stat) => {
      if (!stat.channel_statistics) return;
      const xData = stat.channel_statistics[envelopeXChannel];
      const yData = stat.channel_statistics[envelopeYChannel];
      if (!xData || !yData) return;

      if (!dlcGroups[stat.dlc_number]) {
        dlcGroups[stat.dlc_number] = { x: [], y: [] };
      }

      // Use min, max, and mean as representative points
      dlcGroups[stat.dlc_number].x.push(xData.min, xData.max, xData.mean);
      dlcGroups[stat.dlc_number].y.push(yData.min, yData.max, yData.mean);
    });

    const traces: any[] = [];

    Object.entries(dlcGroups).forEach(([dlc, data]) => {
      traces.push({
        x: data.x,
        y: data.y,
        type: 'scattergl',
        mode: 'markers',
        name: `DLC ${dlc}`,
        marker: {
          color: DLC_COLORS[dlc] ?? '#888',
          size: 5,
          opacity: 0.7,
        },
      });
    });

    // Convex hull overlay: compute from all points
    const allX = Object.values(dlcGroups).flatMap((g) => g.x);
    const allY = Object.values(dlcGroups).flatMap((g) => g.y);

    if (allX.length > 2) {
      const hullPoints = computeConvexHull(
        allX.map((x, i) => [x, allY[i]] as [number, number]),
      );
      if (hullPoints.length > 0) {
        traces.push({
          x: [...hullPoints.map((p) => p[0]), hullPoints[0][0]],
          y: [...hullPoints.map((p) => p[1]), hullPoints[0][1]],
          type: 'scatter',
          mode: 'lines',
          name: 'Envelope',
          line: { color: 'rgba(255,255,255,0.5)', width: 2, dash: 'dash' },
          fill: 'toself',
          fillcolor: 'rgba(255,255,255,0.05)',
        });
      }
    }

    return traces;
  }, [statistics, envelopeXChannel, envelopeYChannel]);

  // ---- Report generation ----
  const handleGenerateReport = useCallback(async () => {
    if (!projectId || !selectedSimId) return;
    setGeneratingReport(true);
    try {
      // Simulate report generation (would call backend)
      await new Promise((resolve) => setTimeout(resolve, 2000));
      toast.success('Report generated successfully');
    } catch {
      toast.error('Failed to generate report');
    } finally {
      setGeneratingReport(false);
    }
  }, [projectId, selectedSimId]);

  if (loading) {
    return (
      <div className="flex h-full items-center justify-center">
        <div className="animate-spin h-8 w-8 rounded-full border-2 border-accent-500 border-t-transparent" />
      </div>
    );
  }

  return (
    <div className="flex h-full flex-col overflow-hidden">
      {/* Simulation selector */}
      <div className="mb-4 flex items-center gap-3">
        <select
          value={selectedSimId ?? ''}
          onChange={(e) => setSelectedSimId(e.target.value)}
          className="input-field max-w-xs"
        >
          {simulations.map((s) => (
            <option key={s.id} value={s.id}>
              {s.name}
            </option>
          ))}
        </select>

        {simulations.length === 0 && (
          <p className="text-sm text-slate-500">No completed simulations available</p>
        )}
      </div>

      {/* Tab navigation */}
      <div className="mb-4 flex border-b border-slate-700">
        {TABS.map((tab) => (
          <button
            key={tab.key}
            onClick={() => setActiveTab(tab.key)}
            className={clsx(
              'flex items-center gap-1.5 border-b-2 px-4 py-2.5 text-sm font-medium transition-colors',
              activeTab === tab.key
                ? 'border-accent-500 text-accent-300'
                : 'border-transparent text-slate-400 hover:text-slate-200',
            )}
          >
            {tab.icon}
            {tab.label}
          </button>
        ))}
      </div>

      {/* Tab content */}
      <div className="flex-1 overflow-auto">
        {activeTab === 'statistics' && (
          <StatisticsTab
            allChannels={allChannels}
            selectedChannel={selectedChannel}
            onSelectChannel={setSelectedChannel}
            channelStats={channelStats}
            aggregatedStats={aggregatedStats}
            boxPlotData={boxPlotData}
          />
        )}

        {activeTab === 'del' && (
          <DELTab
            delResults={delResults}
            delBarData={delBarData}
            mExponent={delMExponent}
            onChangeMExponent={setDelMExponent}
          />
        )}

        {activeTab === 'extreme' && (
          <ExtremeTab extremeEntries={extremeEntries} />
        )}

        {activeTab === 'envelopes' && (
          <EnvelopesTab
            allChannels={allChannels}
            xChannel={envelopeXChannel}
            yChannel={envelopeYChannel}
            onChangeXChannel={setEnvelopeXChannel}
            onChangeYChannel={setEnvelopeYChannel}
            traces={envelopeTraces}
          />
        )}

        {activeTab === 'report' && (
          <ReportTab
            config={reportConfig}
            onChangeConfig={setReportConfig}
            onGenerate={handleGenerateReport}
            generating={generatingReport}
          />
        )}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Statistics Tab
// ---------------------------------------------------------------------------

function StatisticsTab({
  allChannels,
  selectedChannel,
  onSelectChannel,
  channelStats,
  aggregatedStats,
  boxPlotData,
}: {
  allChannels: string[];
  selectedChannel: string;
  onSelectChannel: (ch: string) => void;
  channelStats: Array<{
    dlc: string;
    windSpeed: number;
    min: number;
    max: number;
    mean: number;
    std: number;
    abs_max?: number;
  }>;
  aggregatedStats: { min: number; max: number; mean: number; std: number; abs_max: number } | null;
  boxPlotData: any[];
}) {
  return (
    <div className="space-y-4">
      {/* Channel selector */}
      <div className="flex items-center gap-3">
        <label className="text-sm font-medium text-slate-300">Channel:</label>
        <select
          value={selectedChannel}
          onChange={(e) => onSelectChannel(e.target.value)}
          className="input-field max-w-[200px]"
        >
          {(allChannels.length > 0 ? allChannels : COMMON_CHANNELS).map((ch) => (
            <option key={ch} value={ch}>
              {ch}
            </option>
          ))}
        </select>
      </div>

      {/* Aggregated stats */}
      {aggregatedStats && (
        <div className="grid grid-cols-5 gap-3">
          {[
            { label: 'Min', value: aggregatedStats.min, color: 'text-blue-300' },
            { label: 'Max', value: aggregatedStats.max, color: 'text-red-300' },
            { label: 'Mean', value: aggregatedStats.mean, color: 'text-slate-200' },
            { label: 'Std', value: aggregatedStats.std, color: 'text-amber-300' },
            { label: 'AbsMax', value: aggregatedStats.abs_max, color: 'text-rose-300' },
          ].map((item) => (
            <div
              key={item.label}
              className="rounded-lg border border-slate-700/50 bg-surface-dark p-3 text-center"
            >
              <p className="text-[10px] font-medium uppercase tracking-wider text-slate-500">
                {item.label}
              </p>
              <p className={clsx('mt-1 font-mono text-lg font-semibold', item.color)}>
                {item.value.toFixed(2)}
              </p>
            </div>
          ))}
        </div>
      )}

      {/* Box plot */}
      {boxPlotData.length > 0 && (
        <div className="rounded-xl border border-slate-700/50 bg-surface-dark p-2">
          <Plot
            data={boxPlotData}
            layout={{
              autosize: true,
              height: 300,
              margin: { l: 60, r: 20, t: 30, b: 50 },
              paper_bgcolor: 'rgba(0,0,0,0)',
              plot_bgcolor: 'rgba(15,23,42,0.8)',
              font: { color: '#94a3b8', size: 11 },
              title: {
                text: `${selectedChannel} Distribution by Wind Speed`,
                font: { size: 13, color: '#e2e8f0' },
              },
              xaxis: {
                title: { text: 'Wind Speed', font: { size: 11 } },
                gridcolor: 'rgba(51,65,85,0.4)',
              },
              yaxis: {
                title: { text: selectedChannel, font: { size: 11 } },
                gridcolor: 'rgba(51,65,85,0.4)',
              },
              showlegend: false,
            }}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: '100%' }}
          />
        </div>
      )}

      {/* Statistics table */}
      <div className="rounded-xl border border-slate-700/50 bg-surface-dark-secondary">
        <div className="overflow-auto max-h-[400px]">
          <table className="w-full text-xs">
            <thead className="sticky top-0 bg-surface-dark-secondary z-10">
              <tr className="border-b border-slate-700 text-slate-400">
                <th className="px-3 py-2.5 text-left font-semibold">DLC</th>
                <th className="px-3 py-2.5 text-left font-semibold">Wind Speed</th>
                <th className="px-3 py-2.5 text-right font-semibold">Min</th>
                <th className="px-3 py-2.5 text-right font-semibold">Max</th>
                <th className="px-3 py-2.5 text-right font-semibold">Mean</th>
                <th className="px-3 py-2.5 text-right font-semibold">Std</th>
                <th className="px-3 py-2.5 text-right font-semibold">AbsMax</th>
              </tr>
            </thead>
            <tbody>
              {channelStats.length === 0 && (
                <tr>
                  <td colSpan={7} className="px-4 py-8 text-center text-slate-500">
                    No statistics available
                  </td>
                </tr>
              )}
              {channelStats.map((row, i) => (
                <tr
                  key={i}
                  className="border-b border-slate-800 hover:bg-surface-dark-tertiary/30"
                >
                  <td className="px-3 py-2 font-mono font-medium text-slate-200">
                    {row.dlc}
                  </td>
                  <td className="px-3 py-2 font-mono text-slate-300">
                    {row.windSpeed} m/s
                  </td>
                  <td className="px-3 py-2 text-right font-mono text-blue-300">
                    {row.min.toFixed(2)}
                  </td>
                  <td className="px-3 py-2 text-right font-mono text-red-300">
                    {row.max.toFixed(2)}
                  </td>
                  <td className="px-3 py-2 text-right font-mono text-slate-300">
                    {row.mean.toFixed(2)}
                  </td>
                  <td className="px-3 py-2 text-right font-mono text-slate-400">
                    {row.std.toFixed(2)}
                  </td>
                  <td className="px-3 py-2 text-right font-mono text-rose-300">
                    {(
                      row.abs_max ?? Math.max(Math.abs(row.min), Math.abs(row.max))
                    ).toFixed(2)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// DEL Tab
// ---------------------------------------------------------------------------

function DELTab({
  delResults,
  delBarData,
  mExponent,
  onChangeMExponent,
}: {
  delResults: ResultsDEL | null;
  delBarData: any[];
  mExponent: number;
  onChangeMExponent: (m: number) => void;
}) {
  const delEntries = useMemo(() => {
    if (!delResults?.del_results) return [];
    return Object.entries(delResults.del_results)
      .map(([channel, value]) => ({ channel, value }))
      .sort((a, b) => b.value - a.value);
  }, [delResults]);

  return (
    <div className="space-y-4">
      {/* Woehler exponent selector */}
      <div className="flex items-center gap-4">
        <label className="text-sm font-medium text-slate-300">
          Woehler exponent (m):
        </label>
        <div className="flex gap-2">
          {[
            { m: 4, label: '4 (steel welds)' },
            { m: 10, label: '10 (bolted)' },
            { m: 12, label: '12 (composites)' },
          ].map((opt) => (
            <button
              key={opt.m}
              onClick={() => onChangeMExponent(opt.m)}
              className={clsx(
                'rounded-lg px-3 py-1.5 text-xs font-medium ring-1 ring-inset transition-colors',
                mExponent === opt.m
                  ? 'bg-accent-500/20 text-accent-300 ring-accent-500/40'
                  : 'bg-surface-dark text-slate-400 ring-slate-700 hover:ring-slate-600',
              )}
            >
              {opt.label}
            </button>
          ))}
        </div>
        {delResults && (
          <span className="ml-auto text-xs text-slate-500">
            N_eq = {delResults.n_equivalent.toLocaleString()}
          </span>
        )}
      </div>

      {/* Bar chart */}
      {delBarData.length > 0 && (
        <div className="rounded-xl border border-slate-700/50 bg-surface-dark p-2">
          <Plot
            data={delBarData}
            layout={{
              autosize: true,
              height: 350,
              margin: { l: 60, r: 20, t: 30, b: 100 },
              paper_bgcolor: 'rgba(0,0,0,0)',
              plot_bgcolor: 'rgba(15,23,42,0.8)',
              font: { color: '#94a3b8', size: 11 },
              title: {
                text: 'Damage Equivalent Loads',
                font: { size: 13, color: '#e2e8f0' },
              },
              xaxis: {
                gridcolor: 'rgba(51,65,85,0.4)',
                tickangle: -45,
              },
              yaxis: {
                title: { text: 'DEL', font: { size: 11 } },
                gridcolor: 'rgba(51,65,85,0.4)',
              },
              showlegend: false,
            }}
            config={{ displayModeBar: false, responsive: true }}
            style={{ width: '100%' }}
          />
        </div>
      )}

      {/* Table */}
      <div className="rounded-xl border border-slate-700/50 bg-surface-dark-secondary">
        <div className="overflow-auto max-h-[400px]">
          <table className="w-full text-xs">
            <thead className="sticky top-0 bg-surface-dark-secondary z-10">
              <tr className="border-b border-slate-700 text-slate-400">
                <th className="px-3 py-2.5 text-left font-semibold">Channel</th>
                <th className="px-3 py-2.5 text-right font-semibold">DEL</th>
                <th className="px-3 py-2.5 text-right font-semibold">Units</th>
              </tr>
            </thead>
            <tbody>
              {delEntries.length === 0 && (
                <tr>
                  <td colSpan={3} className="px-4 py-8 text-center text-slate-500">
                    No DEL results available
                  </td>
                </tr>
              )}
              {delEntries.map((entry) => (
                <tr
                  key={entry.channel}
                  className="border-b border-slate-800 hover:bg-surface-dark-tertiary/30"
                >
                  <td className="px-3 py-2 font-mono font-medium text-slate-200">
                    {entry.channel}
                  </td>
                  <td className="px-3 py-2 text-right font-mono text-accent-300">
                    {entry.value.toExponential(3)}
                  </td>
                  <td className="px-3 py-2 text-right text-slate-500">kN/kNm</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Extreme Loads Tab
// ---------------------------------------------------------------------------

function ExtremeTab({
  extremeEntries,
}: {
  extremeEntries: Array<{
    channel: string;
    max: number;
    min: number;
    safety_factor: number;
    design_value: number;
    utilization: number;
  }>;
}) {
  return (
    <div className="space-y-4">
      <div className="rounded-xl border border-slate-700/50 bg-surface-dark-secondary">
        <div className="overflow-auto max-h-[600px]">
          <table className="w-full text-xs">
            <thead className="sticky top-0 bg-surface-dark-secondary z-10">
              <tr className="border-b border-slate-700 text-slate-400">
                <th className="px-3 py-2.5 text-left font-semibold">Channel</th>
                <th className="px-3 py-2.5 text-right font-semibold">Max</th>
                <th className="px-3 py-2.5 text-right font-semibold">Min</th>
                <th className="px-3 py-2.5 text-right font-semibold">Safety Factor</th>
                <th className="px-3 py-2.5 text-right font-semibold">Design Load</th>
                <th className="px-3 py-2.5 text-right font-semibold">Utilization</th>
              </tr>
            </thead>
            <tbody>
              {extremeEntries.length === 0 && (
                <tr>
                  <td colSpan={6} className="px-4 py-8 text-center text-slate-500">
                    No extreme load results available
                  </td>
                </tr>
              )}
              {extremeEntries.map((entry) => {
                const isCritical = entry.utilization > 0.9;
                const isWarning = entry.utilization > 0.7;

                return (
                  <tr
                    key={entry.channel}
                    className={clsx(
                      'border-b border-slate-800 transition-colors',
                      isCritical && 'bg-red-500/10',
                      isWarning && !isCritical && 'bg-amber-500/5',
                      'hover:bg-surface-dark-tertiary/30',
                    )}
                  >
                    <td className="px-3 py-2 font-mono font-medium text-slate-200">
                      <div className="flex items-center gap-2">
                        {entry.channel}
                        {isCritical && (
                          <AlertTriangle size={12} className="text-red-400" />
                        )}
                      </div>
                    </td>
                    <td className="px-3 py-2 text-right font-mono text-red-300">
                      {entry.max.toExponential(3)}
                    </td>
                    <td className="px-3 py-2 text-right font-mono text-blue-300">
                      {entry.min.toExponential(3)}
                    </td>
                    <td className="px-3 py-2 text-right font-mono text-slate-300">
                      {entry.safety_factor.toFixed(2)}
                    </td>
                    <td className="px-3 py-2 text-right font-mono text-amber-300">
                      {entry.design_value.toExponential(3)}
                    </td>
                    <td className="px-3 py-2 text-right">
                      <div className="flex items-center justify-end gap-2">
                        <div className="h-1.5 w-16 overflow-hidden rounded-full bg-slate-700">
                          <div
                            className={clsx(
                              'h-full rounded-full transition-all',
                              isCritical
                                ? 'bg-red-500'
                                : isWarning
                                  ? 'bg-amber-500'
                                  : 'bg-emerald-500',
                            )}
                            style={{
                              width: `${Math.min(100, entry.utilization * 100)}%`,
                            }}
                          />
                        </div>
                        <span
                          className={clsx(
                            'font-mono text-[10px] font-semibold',
                            isCritical
                              ? 'text-red-300'
                              : isWarning
                                ? 'text-amber-300'
                                : 'text-emerald-300',
                          )}
                        >
                          {(entry.utilization * 100).toFixed(1)}%
                        </span>
                      </div>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Load Envelopes Tab
// ---------------------------------------------------------------------------

function EnvelopesTab({
  allChannels,
  xChannel,
  yChannel,
  onChangeXChannel,
  onChangeYChannel,
  traces,
}: {
  allChannels: string[];
  xChannel: string;
  yChannel: string;
  onChangeXChannel: (ch: string) => void;
  onChangeYChannel: (ch: string) => void;
  traces: any[];
}) {
  const channels = allChannels.length > 0 ? allChannels : COMMON_CHANNELS;

  return (
    <div className="space-y-4">
      {/* Axis selectors */}
      <div className="flex items-center gap-6">
        <div className="flex items-center gap-2">
          <label className="text-sm font-medium text-slate-300">X-axis:</label>
          <select
            value={xChannel}
            onChange={(e) => onChangeXChannel(e.target.value)}
            className="input-field max-w-[180px]"
          >
            {channels.map((ch) => (
              <option key={ch} value={ch}>
                {ch}
              </option>
            ))}
          </select>
        </div>
        <div className="flex items-center gap-2">
          <label className="text-sm font-medium text-slate-300">Y-axis:</label>
          <select
            value={yChannel}
            onChange={(e) => onChangeYChannel(e.target.value)}
            className="input-field max-w-[180px]"
          >
            {channels.map((ch) => (
              <option key={ch} value={ch}>
                {ch}
              </option>
            ))}
          </select>
        </div>
      </div>

      {/* Scatter plot */}
      <div className="rounded-xl border border-slate-700/50 bg-surface-dark p-2">
        <Plot
          data={traces}
          layout={{
            autosize: true,
            height: 500,
            margin: { l: 60, r: 20, t: 30, b: 60 },
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(15,23,42,0.8)',
            font: { color: '#94a3b8', size: 11 },
            title: {
              text: `${xChannel} vs ${yChannel} Load Envelope`,
              font: { size: 13, color: '#e2e8f0' },
            },
            xaxis: {
              title: { text: xChannel, font: { size: 11 } },
              gridcolor: 'rgba(51,65,85,0.4)',
              zerolinecolor: 'rgba(51,65,85,0.6)',
            },
            yaxis: {
              title: { text: yChannel, font: { size: 11 } },
              gridcolor: 'rgba(51,65,85,0.4)',
              zerolinecolor: 'rgba(51,65,85,0.6)',
            },
            legend: {
              bgcolor: 'rgba(30,41,59,0.8)',
              bordercolor: 'rgba(51,65,85,0.5)',
              borderwidth: 1,
              font: { size: 10 },
            },
            showlegend: true,
          }}
          config={{ displayModeBar: false, responsive: true }}
          style={{ width: '100%' }}
        />
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Report Tab
// ---------------------------------------------------------------------------

function ReportTab({
  config,
  onChangeConfig,
  onGenerate,
  generating,
}: {
  config: { includeDEL: boolean; includeExtreme: boolean; includeEnvelopes: boolean };
  onChangeConfig: (c: typeof config) => void;
  onGenerate: () => void;
  generating: boolean;
}) {
  return (
    <div className="mx-auto max-w-2xl space-y-6">
      <div className="card">
        <h3 className="mb-4 text-sm font-semibold text-slate-100">Report Configuration</h3>

        <div className="space-y-3">
          {[
            { key: 'includeDEL' as const, label: 'Include Damage Equivalent Loads (DEL)' },
            { key: 'includeExtreme' as const, label: 'Include Extreme Load Analysis' },
            { key: 'includeEnvelopes' as const, label: 'Include Load Envelopes' },
          ].map((item) => (
            <label
              key={item.key}
              className="flex cursor-pointer items-center gap-3 rounded-lg border border-slate-700 p-3 hover:border-slate-600 transition-colors"
            >
              <input
                type="checkbox"
                checked={config[item.key]}
                onChange={(e) =>
                  onChangeConfig({ ...config, [item.key]: e.target.checked })
                }
                className="h-4 w-4 rounded border-slate-500 bg-surface-dark text-accent-500 focus:ring-accent-500 focus:ring-offset-0"
              />
              <span className="text-sm text-slate-200">{item.label}</span>
            </label>
          ))}
        </div>

        <button
          onClick={onGenerate}
          disabled={generating}
          className="btn-primary mt-6 w-full"
        >
          {generating ? (
            <>
              <div className="animate-spin h-4 w-4 rounded-full border-2 border-white border-t-transparent" />
              Generating...
            </>
          ) : (
            <>
              <FileText size={16} />
              Generate Report
            </>
          )}
        </button>
      </div>

      <div className="card">
        <h3 className="mb-4 text-sm font-semibold text-slate-100">Download</h3>
        <div className="flex gap-3">
          <button className="btn-secondary flex-1">
            <Download size={14} />
            HTML Report
          </button>
          <button className="btn-secondary flex-1">
            <Download size={14} />
            PDF Report
          </button>
        </div>
      </div>

      <div className="card">
        <h3 className="mb-4 text-sm font-semibold text-slate-100">Report Preview</h3>
        <div className="flex h-64 items-center justify-center rounded-lg border border-dashed border-slate-600 bg-surface-dark">
          <div className="text-center text-slate-500">
            <FileText size={32} className="mx-auto mb-2 opacity-40" />
            <p className="text-sm">Generate a report to preview</p>
          </div>
        </div>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Convex hull helper (Graham scan)
// ---------------------------------------------------------------------------

function computeConvexHull(points: [number, number][]): [number, number][] {
  if (points.length < 3) return points;

  const pts = [...points].sort((a, b) => a[0] - b[0] || a[1] - b[1]);

  const cross = (o: [number, number], a: [number, number], b: [number, number]) =>
    (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);

  const lower: [number, number][] = [];
  for (const p of pts) {
    while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) {
      lower.pop();
    }
    lower.push(p);
  }

  const upper: [number, number][] = [];
  for (let i = pts.length - 1; i >= 0; i--) {
    const p = pts[i];
    while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) {
      upper.pop();
    }
    upper.push(p);
  }

  lower.pop();
  upper.pop();
  return lower.concat(upper);
}
