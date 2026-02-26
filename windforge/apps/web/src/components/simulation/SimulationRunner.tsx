import { useState, useEffect, useCallback, useRef, useMemo } from 'react';
import { useParams } from 'react-router-dom';
import {
  Play,
  Square,
  RotateCcw,
  Clock,
  Activity,
  AlertTriangle,
  Info,
  XCircle,
  ChevronDown,
  Wifi,
  WifiOff,
  Filter,
} from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';
import Plot from 'react-plotly.js';
import apiClient from '@/api/client';
import StatusBadge from '@/components/common/StatusBadge';
import ProgressBar from '@/components/common/ProgressBar';

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface Simulation {
  id: string;
  project_id: string;
  dlc_definition_id: string;
  turbine_model_id: string;
  name: string;
  status: string;
  total_cases: number;
  completed_cases: number;
  failed_cases: number;
  agent_id: string | null;
  started_at: string | null;
  completed_at: string | null;
  created_at: string;
}

interface SimulationCase {
  id: string;
  simulation_id: string;
  dlc_number: string;
  wind_speed: number;
  seed_number: number;
  yaw_misalignment: number;
  status: string;
  progress_percent: number;
  current_time: number;
  error_message: string | null;
  wall_time_seconds: number | null;
}

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

interface LogEntry {
  time: Date;
  level: 'info' | 'warning' | 'error';
  message: string;
}

interface LiveDataPoint {
  time: number;
  channels: Record<string, number>;
}

// ---------------------------------------------------------------------------
// Available output channels
// ---------------------------------------------------------------------------

const AVAILABLE_CHANNELS = [
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
  'Wind1VelX',
];

const CHANNEL_COLORS: Record<string, string> = {
  GenPwr: '#00b4d8',
  RotSpeed: '#10b981',
  BldPitch1: '#f59e0b',
  TwrBsMxt: '#ef4444',
  TwrBsMyt: '#ec4899',
  RootMxc1: '#8b5cf6',
  RootMyc1: '#6366f1',
  YawBrMxp: '#14b8a6',
  YawBrMyp: '#06b6d4',
  LSShftFxa: '#f97316',
  LSShftTq: '#84cc16',
  TTDspFA: '#a855f7',
  TTDspSS: '#e879f9',
  Wind1VelX: '#22d3ee',
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function formatDuration(seconds: number): string {
  if (seconds < 60) return `${Math.round(seconds)}s`;
  const m = Math.floor(seconds / 60);
  const s = Math.round(seconds % 60);
  if (m < 60) return `${m}m ${s}s`;
  const h = Math.floor(m / 60);
  return `${h}h ${m % 60}m`;
}

function estimateTimeRemaining(sim: Simulation, elapsed: number): string {
  if (sim.completed_cases === 0 || sim.total_cases === 0) return '--';
  const rate = elapsed / sim.completed_cases;
  const remaining = (sim.total_cases - sim.completed_cases) * rate;
  return formatDuration(remaining);
}

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------

export default function SimulationRunner() {
  const { projectId } = useParams<{ projectId: string }>();

  // ---- Simulation state ----
  const [simulations, setSimulations] = useState<Simulation[]>([]);
  const [selectedSimId, setSelectedSimId] = useState<string | null>(null);
  const [simulation, setSimulation] = useState<Simulation | null>(null);
  const [cases, setCases] = useState<SimulationCase[]>([]);
  const [selectedCaseId, setSelectedCaseId] = useState<string | null>(null);
  const [caseStatistics, setCaseStatistics] = useState<ResultsStatistics | null>(null);
  const [loading, setLoading] = useState(true);

  // ---- WebSocket ----
  const wsRef = useRef<WebSocket | null>(null);
  const [wsConnected, setWsConnected] = useState(false);
  const [subscribedCaseId, setSubscribedCaseId] = useState<string | null>(null);

  // ---- Live data ----
  const [liveData, setLiveData] = useState<LiveDataPoint[]>([]);
  const [selectedChannels, setSelectedChannels] = useState<Set<string>>(
    new Set(['GenPwr', 'RotSpeed', 'BldPitch1']),
  );

  // ---- Log ----
  const [logEntries, setLogEntries] = useState<LogEntry[]>([]);
  const [logFilter, setLogFilter] = useState<'all' | 'info' | 'warning' | 'error'>('all');
  const logEndRef = useRef<HTMLDivElement>(null);

  // ---- Timer ----
  const [elapsed, setElapsed] = useState(0);
  const timerRef = useRef<NodeJS.Timeout | null>(null);

  const addLog = useCallback((level: LogEntry['level'], message: string) => {
    setLogEntries((prev) => [...prev.slice(-500), { time: new Date(), level, message }]);
  }, []);

  // ---- Fetch simulations list ----
  useEffect(() => {
    if (!projectId) return;
    setLoading(true);
    apiClient
      .get(`/projects/${projectId}/simulations`)
      .then((res) => {
        const sims = res.data as Simulation[];
        setSimulations(sims);
        if (sims.length > 0 && !selectedSimId) {
          setSelectedSimId(sims[0].id);
        }
      })
      .catch(() => toast.error('Failed to load simulations'))
      .finally(() => setLoading(false));
  }, [projectId]);

  // ---- Fetch selected simulation + cases ----
  useEffect(() => {
    if (!projectId || !selectedSimId) return;

    const fetchSim = () => {
      apiClient
        .get(`/projects/${projectId}/simulations/${selectedSimId}`)
        .then((res) => setSimulation(res.data as Simulation))
        .catch(() => {});
    };

    const fetchCases = () => {
      apiClient
        .get(`/projects/${projectId}/simulations/${selectedSimId}/cases`)
        .then((res) => setCases(res.data as SimulationCase[]))
        .catch(() => {});
    };

    fetchSim();
    fetchCases();

    // Poll while running
    const interval = setInterval(() => {
      if (simulation?.status === 'running') {
        fetchSim();
        fetchCases();
      }
    }, 3000);

    return () => clearInterval(interval);
  }, [projectId, selectedSimId, simulation?.status]);

  // ---- Elapsed timer ----
  useEffect(() => {
    if (simulation?.status === 'running' && simulation.started_at) {
      const startMs = new Date(simulation.started_at).getTime();
      timerRef.current = setInterval(() => {
        setElapsed((Date.now() - startMs) / 1000);
      }, 1000);
    } else {
      if (timerRef.current) clearInterval(timerRef.current);
      if (simulation?.started_at && simulation?.completed_at) {
        setElapsed(
          (new Date(simulation.completed_at).getTime() -
            new Date(simulation.started_at).getTime()) /
            1000,
        );
      }
    }
    return () => {
      if (timerRef.current) clearInterval(timerRef.current);
    };
  }, [simulation?.status, simulation?.started_at, simulation?.completed_at]);

  // ---- WebSocket connection ----
  useEffect(() => {
    if (!selectedSimId || simulation?.status !== 'running') {
      if (wsRef.current) {
        wsRef.current.close();
        wsRef.current = null;
        setWsConnected(false);
      }
      return;
    }

    const ws = new WebSocket(`ws://localhost:8000/ws/${selectedSimId}`);
    wsRef.current = ws;

    ws.onopen = () => {
      setWsConnected(true);
      addLog('info', 'WebSocket connected');
    };

    ws.onclose = () => {
      setWsConnected(false);
      addLog('info', 'WebSocket disconnected');
    };

    ws.onerror = () => {
      addLog('error', 'WebSocket connection error');
    };

    ws.onmessage = (event) => {
      try {
        const msg = JSON.parse(event.data);

        switch (msg.type) {
          case 'case_progress':
            setCases((prev) =>
              prev.map((c) =>
                c.id === msg.case_id
                  ? { ...c, progress_percent: msg.progress_percent, current_time: msg.current_time, status: 'running' }
                  : c,
              ),
            );
            break;

          case 'case_complete':
            setCases((prev) =>
              prev.map((c) =>
                c.id === msg.case_id
                  ? { ...c, status: 'completed', progress_percent: 100 }
                  : c,
              ),
            );
            setSimulation((prev) =>
              prev ? { ...prev, completed_cases: prev.completed_cases + 1 } : prev,
            );
            addLog('info', `Case ${msg.case_id.slice(0, 8)} completed`);
            break;

          case 'case_error':
            setCases((prev) =>
              prev.map((c) =>
                c.id === msg.case_id
                  ? { ...c, status: 'failed', error_message: msg.error }
                  : c,
              ),
            );
            setSimulation((prev) =>
              prev ? { ...prev, failed_cases: prev.failed_cases + 1 } : prev,
            );
            addLog('error', `Case ${msg.case_id.slice(0, 8)} failed: ${msg.error}`);
            break;

          case 'simulation_complete':
            setSimulation((prev) =>
              prev ? { ...prev, status: 'completed' } : prev,
            );
            addLog('info', 'Simulation completed');
            toast.success('Simulation completed');
            break;

          case 'live_data':
            if (msg.case_id === subscribedCaseId) {
              setLiveData((prev) => [
                ...prev.slice(-300),
                { time: msg.time, channels: msg.channels },
              ]);
            }
            break;
        }
      } catch {
        addLog('warning', 'Failed to parse WebSocket message');
      }
    };

    return () => {
      ws.close();
      wsRef.current = null;
      setWsConnected(false);
    };
  }, [selectedSimId, simulation?.status, subscribedCaseId, addLog]);

  // ---- Actions ----
  const handleStart = useCallback(async () => {
    if (!projectId || !selectedSimId) return;
    try {
      await apiClient.post(`/projects/${projectId}/simulations/${selectedSimId}/start`);
      setSimulation((prev) => (prev ? { ...prev, status: 'running' } : prev));
      addLog('info', 'Simulation started');
      toast.success('Simulation started');
    } catch (err: any) {
      toast.error(err?.response?.data?.detail ?? 'Failed to start');
    }
  }, [projectId, selectedSimId, addLog]);

  const handleCancel = useCallback(async () => {
    if (!projectId || !selectedSimId) return;
    try {
      await apiClient.post(`/projects/${projectId}/simulations/${selectedSimId}/cancel`);
      setSimulation((prev) => (prev ? { ...prev, status: 'cancelled' } : prev));
      addLog('warning', 'Simulation cancelled');
      toast('Simulation cancelled');
    } catch (err: any) {
      toast.error(err?.response?.data?.detail ?? 'Failed to cancel');
    }
  }, [projectId, selectedSimId, addLog]);

  const handleRestart = useCallback(async () => {
    if (!projectId) return;
    try {
      const res = await apiClient.post(`/projects/${projectId}/simulations`, {
        dlc_definition_id: simulation?.dlc_definition_id,
        name: `${simulation?.name} (restart)`,
      });
      const newSim = res.data as Simulation;
      setSimulations((prev) => [newSim, ...prev]);
      setSelectedSimId(newSim.id);
      toast.success('New simulation created');
    } catch (err: any) {
      toast.error(err?.response?.data?.detail ?? 'Failed to restart');
    }
  }, [projectId, simulation]);

  // ---- Subscribe to case live data ----
  const handleSubscribe = useCallback(
    (caseId: string) => {
      if (wsRef.current && wsRef.current.readyState === WebSocket.OPEN) {
        wsRef.current.send(JSON.stringify({ command: 'subscribe', case_id: caseId }));
        setSubscribedCaseId(caseId);
        setLiveData([]);
        addLog('info', `Subscribed to live data for case ${caseId.slice(0, 8)}`);
      }
    },
    [addLog],
  );

  // ---- Fetch case statistics when selecting a completed case ----
  useEffect(() => {
    if (!projectId || !selectedSimId || !selectedCaseId) {
      setCaseStatistics(null);
      return;
    }
    const selectedCase = cases.find((c) => c.id === selectedCaseId);
    if (selectedCase?.status !== 'completed') {
      setCaseStatistics(null);
      return;
    }

    apiClient
      .get(`/projects/${projectId}/simulations/${selectedSimId}/results/statistics`)
      .then((res) => {
        const all = res.data as ResultsStatistics[];
        const match = all.find((s) => s.simulation_case_id === selectedCaseId);
        setCaseStatistics(match ?? null);
      })
      .catch(() => {});
  }, [projectId, selectedSimId, selectedCaseId, cases]);

  // ---- Derived ----
  const selectedCase = useMemo(
    () => cases.find((c) => c.id === selectedCaseId) ?? null,
    [cases, selectedCaseId],
  );

  const overallProgress = useMemo(() => {
    if (!simulation || simulation.total_cases === 0) return 0;
    return Math.round(
      ((simulation.completed_cases + simulation.failed_cases) / simulation.total_cases) * 100,
    );
  }, [simulation]);

  const filteredLogs = useMemo(() => {
    if (logFilter === 'all') return logEntries;
    return logEntries.filter((e) => e.level === logFilter);
  }, [logEntries, logFilter]);

  // ---- Plot data ----
  const plotTraces = useMemo(() => {
    if (liveData.length === 0) return [];
    return Array.from(selectedChannels).map((ch) => ({
      x: liveData.map((d) => d.time),
      y: liveData.map((d) => d.channels[ch] ?? 0),
      type: 'scattergl' as const,
      mode: 'lines' as const,
      name: ch,
      line: { color: CHANNEL_COLORS[ch] ?? '#888', width: 1.5 },
    }));
  }, [liveData, selectedChannels]);

  // ---- Sort cases for display ----
  const sortedCases = useMemo(() => {
    return [...cases].sort((a, b) => {
      if (a.dlc_number !== b.dlc_number) return a.dlc_number.localeCompare(b.dlc_number);
      if (a.wind_speed !== b.wind_speed) return a.wind_speed - b.wind_speed;
      return a.seed_number - b.seed_number;
    });
  }, [cases]);

  // Scroll log to bottom
  useEffect(() => {
    logEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [logEntries.length]);

  if (loading) {
    return (
      <div className="flex h-full items-center justify-center">
        <div className="animate-spin h-8 w-8 rounded-full border-2 border-accent-500 border-t-transparent" />
      </div>
    );
  }

  return (
    <div className="flex h-full flex-col overflow-hidden">
      {/* ====== Top Bar ====== */}
      <div className="mb-4 flex flex-wrap items-center gap-3">
        {/* Simulation selector */}
        <select
          value={selectedSimId ?? ''}
          onChange={(e) => {
            setSelectedSimId(e.target.value);
            setSelectedCaseId(null);
            setLiveData([]);
          }}
          className="input-field max-w-xs"
        >
          {simulations.map((s) => (
            <option key={s.id} value={s.id}>
              {s.name}
            </option>
          ))}
        </select>

        {simulation && <StatusBadge status={simulation.status} />}

        {/* Connection indicator */}
        <span
          className={clsx(
            'flex items-center gap-1 text-xs',
            wsConnected ? 'text-emerald-400' : 'text-slate-500',
          )}
        >
          {wsConnected ? <Wifi size={12} /> : <WifiOff size={12} />}
          {wsConnected ? 'Live' : 'Offline'}
        </span>

        {/* Actions */}
        <div className="ml-auto flex gap-2">
          {simulation?.status === 'pending' && (
            <button onClick={handleStart} className="btn-primary text-xs">
              <Play size={14} /> Start
            </button>
          )}
          {simulation?.status === 'running' && (
            <button onClick={handleCancel} className="btn-danger text-xs">
              <Square size={14} /> Cancel
            </button>
          )}
          {(simulation?.status === 'completed' ||
            simulation?.status === 'failed' ||
            simulation?.status === 'cancelled') && (
            <button onClick={handleRestart} className="btn-secondary text-xs">
              <RotateCcw size={14} /> Restart
            </button>
          )}
        </div>
      </div>

      {/* Progress overview */}
      {simulation && (
        <div className="mb-4 rounded-xl border border-slate-700/50 bg-surface-dark-secondary p-4">
          <div className="mb-2 flex items-center justify-between text-sm">
            <span className="text-slate-400">
              {simulation.completed_cases} / {simulation.total_cases} cases completed
              {simulation.failed_cases > 0 && (
                <span className="ml-2 text-red-400">
                  ({simulation.failed_cases} failed)
                </span>
              )}
            </span>
            <div className="flex items-center gap-4 text-xs text-slate-400">
              <span className="flex items-center gap-1">
                <Clock size={12} />
                Elapsed: {formatDuration(elapsed)}
              </span>
              {simulation.status === 'running' && (
                <span className="flex items-center gap-1">
                  <Activity size={12} />
                  ETA: {estimateTimeRemaining(simulation, elapsed)}
                </span>
              )}
            </div>
          </div>
          <ProgressBar
            value={overallProgress}
            size="md"
            color={
              simulation.status === 'failed'
                ? 'red'
                : simulation.status === 'completed'
                  ? 'green'
                  : 'accent'
            }
            showLabel
          />
        </div>
      )}

      {/* ====== Main Content ====== */}
      <div className="flex flex-1 gap-4 overflow-hidden">
        {/* Case Matrix Grid */}
        <div className="flex-1 overflow-auto rounded-xl border border-slate-700/50 bg-surface-dark-secondary">
          <table className="w-full text-xs">
            <thead className="sticky top-0 z-10 bg-surface-dark-secondary">
              <tr className="border-b border-slate-700 text-slate-400">
                <th className="px-3 py-2.5 text-left font-semibold">DLC</th>
                <th className="px-3 py-2.5 text-left font-semibold">Wind (m/s)</th>
                <th className="px-3 py-2.5 text-left font-semibold">Seed</th>
                <th className="px-3 py-2.5 text-left font-semibold">Yaw</th>
                <th className="px-3 py-2.5 text-left font-semibold">Status</th>
                <th className="px-3 py-2.5 text-left font-semibold w-40">Progress</th>
                <th className="px-3 py-2.5 text-right font-semibold">Time</th>
              </tr>
            </thead>
            <tbody>
              {sortedCases.length === 0 && (
                <tr>
                  <td colSpan={7} className="px-4 py-12 text-center text-slate-500">
                    No cases available
                  </td>
                </tr>
              )}
              {sortedCases.map((c) => (
                <tr
                  key={c.id}
                  onClick={() => setSelectedCaseId(c.id)}
                  className={clsx(
                    'cursor-pointer border-b border-slate-800 transition-colors',
                    selectedCaseId === c.id
                      ? 'bg-accent-500/10'
                      : 'hover:bg-surface-dark-tertiary/30',
                  )}
                >
                  <td className="px-3 py-2 font-mono font-medium text-slate-200">
                    {c.dlc_number}
                  </td>
                  <td className="px-3 py-2 font-mono text-slate-300">{c.wind_speed}</td>
                  <td className="px-3 py-2 font-mono text-slate-300">{c.seed_number}</td>
                  <td className="px-3 py-2 font-mono text-slate-300">
                    {c.yaw_misalignment > 0 ? `+${c.yaw_misalignment}` : c.yaw_misalignment}°
                  </td>
                  <td className="px-3 py-2">
                    <StatusBadge status={c.status} size="sm" />
                  </td>
                  <td className="px-3 py-2">
                    {c.status === 'running' ? (
                      <ProgressBar value={c.progress_percent} size="sm" color="amber" animated />
                    ) : c.status === 'completed' ? (
                      <ProgressBar value={100} size="sm" color="green" animated={false} />
                    ) : (
                      <ProgressBar value={0} size="sm" color="accent" animated={false} />
                    )}
                  </td>
                  <td className="px-3 py-2 text-right font-mono text-slate-400">
                    {c.wall_time_seconds != null ? formatDuration(c.wall_time_seconds) : '--'}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>

        {/* ====== Right Panel - Live Case View ====== */}
        <div className="w-96 shrink-0 overflow-y-auto rounded-xl border border-slate-700/50 bg-surface-dark-secondary p-4">
          {!selectedCase && (
            <div className="flex h-full flex-col items-center justify-center text-slate-500">
              <Activity size={32} className="mb-3 opacity-40" />
              <p className="text-sm">Select a case to view details</p>
            </div>
          )}

          {selectedCase && selectedCase.status === 'running' && (
            <>
              <div className="mb-4">
                <h3 className="mb-1 text-sm font-semibold text-slate-100">
                  Live View: DLC {selectedCase.dlc_number}
                </h3>
                <p className="text-xs text-slate-400">
                  Wind: {selectedCase.wind_speed} m/s | Seed: {selectedCase.seed_number} | Yaw:{' '}
                  {selectedCase.yaw_misalignment}°
                </p>
                <div className="mt-2 flex items-center gap-2">
                  <ProgressBar value={selectedCase.progress_percent} size="sm" color="amber" />
                  <span className="whitespace-nowrap text-xs font-mono text-amber-300">
                    {selectedCase.progress_percent.toFixed(0)}%
                  </span>
                </div>
                <p className="mt-1 text-xs text-slate-500">
                  Sim time: {selectedCase.current_time.toFixed(1)}s
                </p>
              </div>

              {/* Subscribe button */}
              {subscribedCaseId !== selectedCase.id && (
                <button
                  onClick={() => handleSubscribe(selectedCase.id)}
                  disabled={!wsConnected}
                  className="btn-secondary mb-4 w-full text-xs"
                >
                  <Activity size={14} />
                  Subscribe to Live Data
                </button>
              )}

              {/* Channel selector */}
              <div className="mb-3">
                <p className="mb-1.5 text-xs font-semibold text-slate-400 uppercase tracking-wider">
                  Channels
                </p>
                <div className="flex flex-wrap gap-1">
                  {AVAILABLE_CHANNELS.map((ch) => (
                    <button
                      key={ch}
                      onClick={() =>
                        setSelectedChannels((prev) => {
                          const next = new Set(prev);
                          if (next.has(ch)) {
                            if (next.size > 1) next.delete(ch);
                          } else {
                            next.add(ch);
                          }
                          return next;
                        })
                      }
                      className={clsx(
                        'rounded px-2 py-0.5 text-[10px] font-mono font-medium ring-1 ring-inset transition-colors',
                        selectedChannels.has(ch)
                          ? 'bg-accent-500/20 text-accent-300 ring-accent-500/40'
                          : 'bg-surface-dark text-slate-500 ring-slate-700 hover:ring-slate-600',
                      )}
                    >
                      {ch}
                    </button>
                  ))}
                </div>
              </div>

              {/* Live plot */}
              {subscribedCaseId === selectedCase.id && liveData.length > 0 && (
                <div className="rounded-lg border border-slate-700 bg-surface-dark p-1">
                  <Plot
                    data={plotTraces}
                    layout={{
                      autosize: true,
                      height: 280,
                      margin: { l: 50, r: 10, t: 10, b: 40 },
                      paper_bgcolor: 'rgba(0,0,0,0)',
                      plot_bgcolor: 'rgba(15,23,42,0.8)',
                      font: { color: '#94a3b8', size: 10 },
                      xaxis: {
                        title: { text: 'Time (s)', font: { size: 10 } },
                        gridcolor: 'rgba(51,65,85,0.4)',
                        zerolinecolor: 'rgba(51,65,85,0.4)',
                      },
                      yaxis: {
                        gridcolor: 'rgba(51,65,85,0.4)',
                        zerolinecolor: 'rgba(51,65,85,0.4)',
                      },
                      legend: { orientation: 'h', y: -0.25, font: { size: 9 } },
                      showlegend: true,
                    }}
                    config={{ displayModeBar: false, responsive: true }}
                    style={{ width: '100%' }}
                  />
                </div>
              )}

              {subscribedCaseId === selectedCase.id && liveData.length === 0 && (
                <p className="py-8 text-center text-xs text-slate-500">
                  Waiting for live data...
                </p>
              )}
            </>
          )}

          {selectedCase && selectedCase.status === 'completed' && (
            <>
              <div className="mb-4">
                <h3 className="mb-1 text-sm font-semibold text-slate-100">
                  Completed: DLC {selectedCase.dlc_number}
                </h3>
                <p className="text-xs text-slate-400">
                  Wind: {selectedCase.wind_speed} m/s | Seed: {selectedCase.seed_number} | Yaw:{' '}
                  {selectedCase.yaw_misalignment}°
                </p>
                {selectedCase.wall_time_seconds != null && (
                  <p className="mt-1 text-xs text-slate-500">
                    Wall time: {formatDuration(selectedCase.wall_time_seconds)}
                  </p>
                )}
              </div>

              {/* Statistics table */}
              {caseStatistics?.channel_statistics && (
                <div className="overflow-auto rounded-lg border border-slate-700">
                  <table className="w-full text-[10px]">
                    <thead>
                      <tr className="border-b border-slate-700 bg-surface-dark text-slate-400">
                        <th className="px-2 py-1.5 text-left font-semibold">Channel</th>
                        <th className="px-2 py-1.5 text-right font-semibold">Min</th>
                        <th className="px-2 py-1.5 text-right font-semibold">Max</th>
                        <th className="px-2 py-1.5 text-right font-semibold">Mean</th>
                        <th className="px-2 py-1.5 text-right font-semibold">Std</th>
                      </tr>
                    </thead>
                    <tbody>
                      {Object.entries(caseStatistics.channel_statistics).map(
                        ([channel, stats]) => (
                          <tr
                            key={channel}
                            className="border-b border-slate-800 hover:bg-surface-dark-tertiary/30"
                          >
                            <td className="px-2 py-1 font-mono text-slate-200">{channel}</td>
                            <td className="px-2 py-1 text-right font-mono text-blue-300">
                              {stats.min.toFixed(2)}
                            </td>
                            <td className="px-2 py-1 text-right font-mono text-red-300">
                              {stats.max.toFixed(2)}
                            </td>
                            <td className="px-2 py-1 text-right font-mono text-slate-300">
                              {stats.mean.toFixed(2)}
                            </td>
                            <td className="px-2 py-1 text-right font-mono text-slate-400">
                              {stats.std.toFixed(2)}
                            </td>
                          </tr>
                        ),
                      )}
                    </tbody>
                  </table>
                </div>
              )}

              {!caseStatistics && (
                <p className="py-8 text-center text-xs text-slate-500">
                  No statistics available for this case
                </p>
              )}
            </>
          )}

          {selectedCase && selectedCase.status === 'failed' && (
            <div className="rounded-lg border border-red-500/30 bg-red-500/10 p-4">
              <div className="mb-2 flex items-center gap-2 text-sm font-semibold text-red-300">
                <XCircle size={16} />
                Case Failed
              </div>
              <p className="text-xs text-red-200/80">
                DLC {selectedCase.dlc_number} | Wind: {selectedCase.wind_speed} m/s
              </p>
              {selectedCase.error_message && (
                <pre className="mt-3 overflow-auto rounded bg-red-900/30 p-2 font-mono text-[10px] text-red-200">
                  {selectedCase.error_message}
                </pre>
              )}
            </div>
          )}

          {selectedCase &&
            selectedCase.status !== 'running' &&
            selectedCase.status !== 'completed' &&
            selectedCase.status !== 'failed' && (
              <div className="flex h-32 items-center justify-center text-sm text-slate-500">
                Case is {selectedCase.status}
              </div>
            )}
        </div>
      </div>

      {/* ====== Bottom Log Panel ====== */}
      <div className="mt-4 flex h-40 flex-col rounded-xl border border-slate-700/50 bg-surface-dark-secondary">
        {/* Log header */}
        <div className="flex items-center gap-2 border-b border-slate-700 px-4 py-2">
          <span className="text-xs font-semibold text-slate-300">Console</span>
          <div className="ml-auto flex gap-1">
            {(['all', 'info', 'warning', 'error'] as const).map((level) => (
              <button
                key={level}
                onClick={() => setLogFilter(level)}
                className={clsx(
                  'rounded px-2 py-0.5 text-[10px] font-medium transition-colors',
                  logFilter === level
                    ? 'bg-accent-500/20 text-accent-300'
                    : 'text-slate-500 hover:text-slate-300',
                )}
              >
                {level === 'all' && <Filter size={10} className="inline mr-0.5" />}
                {level}
              </button>
            ))}
          </div>
        </div>

        {/* Log entries */}
        <div className="flex-1 overflow-y-auto px-4 py-1 font-mono text-[11px]">
          {filteredLogs.length === 0 && (
            <p className="py-4 text-center text-slate-600">No log entries</p>
          )}
          {filteredLogs.map((entry, i) => (
            <div key={i} className="flex gap-2 py-0.5">
              <span className="shrink-0 text-slate-600">
                {entry.time.toLocaleTimeString()}
              </span>
              <span
                className={clsx(
                  'shrink-0 w-12',
                  entry.level === 'info' && 'text-blue-400',
                  entry.level === 'warning' && 'text-amber-400',
                  entry.level === 'error' && 'text-red-400',
                )}
              >
                {entry.level === 'info' && <Info size={11} className="inline" />}
                {entry.level === 'warning' && <AlertTriangle size={11} className="inline" />}
                {entry.level === 'error' && <XCircle size={11} className="inline" />}
                {' ' + entry.level.toUpperCase()}
              </span>
              <span className="text-slate-300">{entry.message}</span>
            </div>
          ))}
          <div ref={logEndRef} />
        </div>
      </div>
    </div>
  );
}
