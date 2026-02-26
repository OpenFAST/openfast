import React, { useCallback, useEffect, useMemo, useState } from 'react';
import { useParams } from 'react-router-dom';
import axios from 'axios';
import {
  Save,
  RotateCcw,
  Plus,
  ChevronDown,
  Loader2,
  AlertCircle,
  Cpu,
  Settings,
} from 'lucide-react';
import clsx from 'clsx';
import FilePreview from '@/components/common/FilePreview';
import PlotPanel from '@/components/common/PlotPanel';

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface Controller {
  id: string;
  project_id: string;
  name: string;
  version: number;
  controller_type: string;
  pcmode: number;
  vscontrl: number;
  parameters: Record<string, any> | null;
  dll_filename: string | null;
  dll_procname: string | null;
  is_active: boolean;
  created_at: string;
}

const API_BASE = '/api/v1';

const CONTROLLER_TYPES = [
  { value: 'baseline', label: 'Baseline PI' },
  { value: 'rosco', label: 'ROSCO' },
  { value: 'discon', label: 'DISCON DLL' },
  { value: 'bladed', label: 'Bladed DLL' },
  { value: 'custom', label: 'Custom' },
];

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function ParamInput({
  label,
  value,
  unit,
  onChange,
  type = 'number',
}: {
  label: string;
  value: any;
  unit?: string;
  onChange: (v: any) => void;
  type?: 'number' | 'text';
}) {
  return (
    <div className="flex items-center gap-2">
      <label className="text-xs text-gray-400 w-36 shrink-0 truncate" title={label}>
        {label}
      </label>
      <input
        type={type}
        step={type === 'number' ? 'any' : undefined}
        value={value ?? (type === 'number' ? 0 : '')}
        onChange={(e) =>
          onChange(type === 'number' ? parseFloat(e.target.value) || 0 : e.target.value)
        }
        className="flex-1 px-2 py-1 bg-gray-800 border border-gray-700 rounded
                   text-xs font-mono text-gray-200 focus:border-blue-500
                   focus:outline-none focus:ring-1 focus:ring-blue-500/50
                   [appearance:textfield]
                   [&::-webkit-outer-spin-button]:appearance-none
                   [&::-webkit-inner-spin-button]:appearance-none"
      />
      {unit && <span className="text-xs text-gray-500 w-14 shrink-0">{unit}</span>}
    </div>
  );
}

function SectionHeader({ title, icon: Icon }: { title: string; icon?: any }) {
  return (
    <div className="flex items-center gap-2 mb-2 mt-4 first:mt-0">
      {Icon && <Icon size={13} className="text-gray-500" />}
      <span className="text-xs font-bold text-gray-300 uppercase tracking-wider">{title}</span>
      <div className="flex-1 border-b border-gray-800" />
    </div>
  );
}

// ---------------------------------------------------------------------------
// Gain schedule table
// ---------------------------------------------------------------------------

function GainScheduleTable({
  schedule,
  onChange,
}: {
  schedule: { wind_speed: number; kp: number; ki: number }[];
  onChange: (s: { wind_speed: number; kp: number; ki: number }[]) => void;
}) {
  return (
    <div className="border border-gray-800 rounded overflow-auto max-h-[200px]">
      <table className="w-full text-xs">
        <thead className="bg-gray-800 sticky top-0">
          <tr>
            <th className="px-2 py-1.5 text-left text-gray-400 font-semibold">Wind Speed [m/s]</th>
            <th className="px-2 py-1.5 text-left text-gray-400 font-semibold">KP [s]</th>
            <th className="px-2 py-1.5 text-left text-gray-400 font-semibold">KI [-]</th>
            <th className="w-8" />
          </tr>
        </thead>
        <tbody>
          {schedule.map((row, i) => (
            <tr key={i} className={i % 2 === 0 ? 'bg-gray-900' : 'bg-gray-900/60'}>
              <td className="px-1 py-0.5">
                <input
                  type="number"
                  step="any"
                  value={row.wind_speed}
                  onChange={(e) => {
                    const next = [...schedule];
                    next[i] = { ...next[i], wind_speed: parseFloat(e.target.value) || 0 };
                    onChange(next);
                  }}
                  className="w-full px-1 py-0.5 bg-gray-800 border border-gray-700 rounded
                             text-xs font-mono text-gray-200
                             [appearance:textfield]
                             [&::-webkit-outer-spin-button]:appearance-none
                             [&::-webkit-inner-spin-button]:appearance-none"
                />
              </td>
              <td className="px-1 py-0.5">
                <input
                  type="number"
                  step="any"
                  value={row.kp}
                  onChange={(e) => {
                    const next = [...schedule];
                    next[i] = { ...next[i], kp: parseFloat(e.target.value) || 0 };
                    onChange(next);
                  }}
                  className="w-full px-1 py-0.5 bg-gray-800 border border-gray-700 rounded
                             text-xs font-mono text-gray-200
                             [appearance:textfield]
                             [&::-webkit-outer-spin-button]:appearance-none
                             [&::-webkit-inner-spin-button]:appearance-none"
                />
              </td>
              <td className="px-1 py-0.5">
                <input
                  type="number"
                  step="any"
                  value={row.ki}
                  onChange={(e) => {
                    const next = [...schedule];
                    next[i] = { ...next[i], ki: parseFloat(e.target.value) || 0 };
                    onChange(next);
                  }}
                  className="w-full px-1 py-0.5 bg-gray-800 border border-gray-700 rounded
                             text-xs font-mono text-gray-200
                             [appearance:textfield]
                             [&::-webkit-outer-spin-button]:appearance-none
                             [&::-webkit-inner-spin-button]:appearance-none"
                />
              </td>
              <td className="px-1 py-0.5 text-center">
                <button
                  onClick={() => onChange(schedule.filter((_, j) => j !== i))}
                  className="text-gray-600 hover:text-red-400 transition-colors"
                  title="Remove row"
                >
                  &times;
                </button>
              </td>
            </tr>
          ))}
        </tbody>
      </table>
      <button
        onClick={() => {
          const lastWS = schedule.length > 0 ? schedule[schedule.length - 1].wind_speed + 1 : 5;
          onChange([...schedule, { wind_speed: lastWS, kp: 0, ki: 0 }]);
        }}
        className="w-full py-1 text-xs text-gray-500 hover:text-blue-400 hover:bg-gray-800
                   transition-colors border-t border-gray-800"
      >
        + Add Row
      </button>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Torque-speed curve visualization
// ---------------------------------------------------------------------------

function TorqueSpeedPlot({ controller }: { controller: Controller }) {
  const params = controller.parameters ?? {};
  const ratedGenSpeed = params.rated_gen_speed ?? 1173.7;
  const rgn2K = params.rgn2k ?? 0.0255764;
  const ratedTorque = params.rated_torque ?? 43093.55;
  const minGenSpeed = params.min_gen_speed ?? 670;
  const maxGenSpeed = params.max_gen_speed ?? 1173.7;
  const rgn15StartSpeed = minGenSpeed * 1.0;
  const rgn2StartSpeed = minGenSpeed * 1.1;

  // Build curve points
  const speeds: number[] = [];
  const torques: number[] = [];
  const colors: string[] = [];
  const regionLabels: { speed: number; torque: number; label: string }[] = [];

  // Region 1 (below cut-in)
  for (let s = 0; s < rgn15StartSpeed; s += 20) {
    speeds.push(s);
    torques.push(0);
    colors.push('#6b7280');
  }
  regionLabels.push({ speed: rgn15StartSpeed / 2, torque: 0, label: 'Rgn 1' });

  // Region 1.5 (linear transition)
  const rgn15Steps = 10;
  for (let i = 0; i <= rgn15Steps; i++) {
    const s = rgn15StartSpeed + (i / rgn15Steps) * (rgn2StartSpeed - rgn15StartSpeed);
    const t = rgn2K * s * s * (i / rgn15Steps);
    speeds.push(s);
    torques.push(t);
    colors.push('#f59e0b');
  }
  regionLabels.push({
    speed: (rgn15StartSpeed + rgn2StartSpeed) / 2,
    torque: rgn2K * rgn2StartSpeed * rgn2StartSpeed * 0.3,
    label: 'Rgn 1.5',
  });

  // Region 2 (quadratic)
  const rgn25Start = ratedGenSpeed * 0.95;
  for (let s = rgn2StartSpeed; s <= rgn25Start; s += 10) {
    speeds.push(s);
    torques.push(rgn2K * s * s);
    colors.push('#3b82f6');
  }
  regionLabels.push({
    speed: (rgn2StartSpeed + rgn25Start) / 2,
    torque: rgn2K * ((rgn2StartSpeed + rgn25Start) / 2) ** 2,
    label: 'Rgn 2',
  });

  // Region 2.5 (transition to rated)
  const rgn25Steps = 10;
  for (let i = 0; i <= rgn25Steps; i++) {
    const s = rgn25Start + (i / rgn25Steps) * (ratedGenSpeed - rgn25Start);
    const t2 = rgn2K * s * s;
    const t = t2 + (ratedTorque - t2) * (i / rgn25Steps);
    speeds.push(s);
    torques.push(t);
    colors.push('#8b5cf6');
  }
  regionLabels.push({
    speed: (rgn25Start + ratedGenSpeed) / 2,
    torque: ratedTorque * 0.85,
    label: 'Rgn 2.5',
  });

  // Region 3 (constant torque, pitch control)
  for (let s = ratedGenSpeed; s <= maxGenSpeed * 1.05; s += 10) {
    speeds.push(s);
    torques.push(ratedTorque);
    colors.push('#ef4444');
  }
  regionLabels.push({
    speed: ratedGenSpeed + (maxGenSpeed * 1.05 - ratedGenSpeed) / 2,
    torque: ratedTorque,
    label: 'Rgn 3',
  });

  const plotData = [
    {
      x: speeds,
      y: torques,
      type: 'scatter' as const,
      mode: 'lines' as const,
      line: { color: '#3b82f6', width: 2 },
      name: 'Torque-Speed',
    },
  ];

  // Add region annotations
  const annotations = regionLabels.map((rl) => ({
    x: rl.speed,
    y: rl.torque,
    text: rl.label,
    showarrow: false,
    font: { size: 9, color: '#9ca3af' },
    yshift: 15,
  }));

  return (
    <PlotPanel
      data={plotData}
      layout={{
        xaxis: { title: { text: 'Gen Speed (rpm)' } },
        yaxis: { title: { text: 'Gen Torque (N-m)' } },
        annotations,
      }}
      title="Operating Regions"
      height={300}
    />
  );
}

// ---------------------------------------------------------------------------
// Main ControllerDesigner
// ---------------------------------------------------------------------------

export default function ControllerDesigner() {
  const { projectId } = useParams<{ projectId: string }>();

  const [controllers, setControllers] = useState<Controller[]>([]);
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [controller, setController] = useState<Controller | null>(null);
  const [originalController, setOriginalController] = useState<Controller | null>(null);
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [rightTab, setRightTab] = useState<'regions' | 'pitch_schedule' | 'preview'>('regions');
  const [filePreview, setFilePreview] = useState<string | null>(null);
  const [previewLoading, setPreviewLoading] = useState(false);
  const [dropdownOpen, setDropdownOpen] = useState(false);

  const params = controller?.parameters ?? {};

  const isDirty = useMemo(
    () => JSON.stringify(controller) !== JSON.stringify(originalController),
    [controller, originalController]
  );

  // ---------------------------------------------------------------------------
  // Fetching
  // ---------------------------------------------------------------------------

  const fetchControllers = useCallback(async () => {
    if (!projectId) return;
    try {
      setLoading(true);
      const resp = await axios.get<Controller[]>(
        `${API_BASE}/projects/${projectId}/controllers`
      );
      setControllers(resp.data);
      if (resp.data.length > 0 && !selectedId) {
        setSelectedId(resp.data[0].id);
      }
    } catch (err: any) {
      setError(err.message || 'Failed to load controllers');
    } finally {
      setLoading(false);
    }
  }, [projectId, selectedId]);

  const fetchController = useCallback(
    async (controllerId: string) => {
      if (!projectId) return;
      try {
        setLoading(true);
        const resp = await axios.get<Controller>(
          `${API_BASE}/projects/${projectId}/controllers/${controllerId}`
        );
        setController(resp.data);
        setOriginalController(JSON.parse(JSON.stringify(resp.data)));
        setError(null);
      } catch (err: any) {
        setError(err.message || 'Failed to load controller');
      } finally {
        setLoading(false);
      }
    },
    [projectId]
  );

  const fetchPreview = useCallback(async () => {
    if (!projectId || !selectedId) return;
    try {
      setPreviewLoading(true);
      const resp = await axios.get<string>(
        `${API_BASE}/projects/${projectId}/controllers/${selectedId}/preview`
      );
      setFilePreview(resp.data);
    } catch {
      setFilePreview(null);
    } finally {
      setPreviewLoading(false);
    }
  }, [projectId, selectedId]);

  useEffect(() => {
    fetchControllers();
  }, [fetchControllers]);

  useEffect(() => {
    if (selectedId) fetchController(selectedId);
  }, [selectedId, fetchController]);

  useEffect(() => {
    if (rightTab === 'preview') fetchPreview();
  }, [rightTab, fetchPreview]);

  // ---------------------------------------------------------------------------
  // Actions
  // ---------------------------------------------------------------------------

  const handleSave = useCallback(async () => {
    if (!projectId || !controller) return;
    try {
      setSaving(true);
      if (controller.id) {
        await axios.put(
          `${API_BASE}/projects/${projectId}/controllers/${controller.id}`,
          controller
        );
      } else {
        const resp = await axios.post<Controller>(
          `${API_BASE}/projects/${projectId}/controllers`,
          controller
        );
        setSelectedId(resp.data.id);
      }
      setOriginalController(JSON.parse(JSON.stringify(controller)));
      await fetchControllers();
      setError(null);
    } catch (err: any) {
      setError(err.message || 'Failed to save controller');
    } finally {
      setSaving(false);
    }
  }, [projectId, controller, fetchControllers]);

  const handleRevert = useCallback(() => {
    if (originalController) setController(JSON.parse(JSON.stringify(originalController)));
  }, [originalController]);

  const handleNewController = useCallback(() => {
    const newCtrl: Controller = {
      id: '',
      project_id: projectId || '',
      name: 'New Controller',
      version: 1,
      controller_type: 'baseline',
      pcmode: 0,
      vscontrl: 1,
      parameters: {
        tpcon: 0,
        pitch_min: 0,
        pitch_max: 90,
        pitch_rate_limit: 8,
        kp_pitch: 0.01882681,
        ki_pitch: 0.008068634,
        rated_gen_speed: 1173.7,
        rgn2k: 0.0255764,
        rated_torque: 43093.55,
        min_gen_torque: 0,
        max_gen_torque: 47402.91,
        min_gen_speed: 670,
        max_gen_speed: 1173.7,
        gen_model: 1,
        gen_eff: 94.4,
        gain_schedule: [
          { wind_speed: 11.4, kp: 0.01882681, ki: 0.008068634 },
          { wind_speed: 13, kp: 0.01350, ki: 0.005790 },
          { wind_speed: 15, kp: 0.00980, ki: 0.004200 },
          { wind_speed: 18, kp: 0.00670, ki: 0.002870 },
          { wind_speed: 21, kp: 0.00490, ki: 0.002100 },
          { wind_speed: 25, kp: 0.00370, ki: 0.001590 },
        ],
      },
      dll_filename: null,
      dll_procname: null,
      is_active: true,
      created_at: new Date().toISOString(),
    };
    setController(newCtrl);
    setOriginalController(null);
    setSelectedId(null);
  }, [projectId]);

  const updateField = useCallback(
    <K extends keyof Controller>(key: K, value: Controller[K]) => {
      if (!controller) return;
      setController({ ...controller, [key]: value });
    },
    [controller]
  );

  const updateParam = useCallback(
    (key: string, value: any) => {
      if (!controller) return;
      setController({
        ...controller,
        parameters: { ...controller.parameters, [key]: value },
      });
    },
    [controller]
  );

  // ---------------------------------------------------------------------------
  // Plot data for gain schedule
  // ---------------------------------------------------------------------------

  const gainSchedulePlot = useMemo(() => {
    const schedule = (params.gain_schedule as { wind_speed: number; kp: number; ki: number }[]) ?? [];
    const sorted = [...schedule].sort((a, b) => a.wind_speed - b.wind_speed);

    return {
      data: [
        {
          x: sorted.map((s) => s.wind_speed),
          y: sorted.map((s) => s.kp),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'KP',
          line: { color: '#3b82f6', width: 2 },
          marker: { size: 5 },
          yaxis: 'y' as const,
        },
        {
          x: sorted.map((s) => s.wind_speed),
          y: sorted.map((s) => s.ki),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'KI',
          line: { color: '#ef4444', width: 2, dash: 'dash' as const },
          marker: { size: 5 },
          yaxis: 'y2' as const,
        },
      ],
      layout: {
        xaxis: { title: { text: 'Wind Speed (m/s)' } },
        yaxis: { title: { text: 'KP (s)' }, side: 'left' as const },
        yaxis2: {
          title: { text: 'KI (-)' },
          side: 'right' as const,
          overlaying: 'y' as const,
        },
      },
    };
  }, [params.gain_schedule]);

  // ---------------------------------------------------------------------------
  // Render
  // ---------------------------------------------------------------------------

  if (loading && !controller) {
    return (
      <div className="flex items-center justify-center h-full bg-gray-950 text-gray-400">
        <Loader2 className="animate-spin mr-2" size={20} />
        Loading controller data...
      </div>
    );
  }

  return (
    <div className="flex flex-col h-full bg-gray-950 text-gray-200">
      {/* Top bar */}
      <div className="flex items-center gap-3 px-4 py-2 bg-gray-900 border-b border-gray-800 shrink-0">
        <Cpu size={16} className="text-purple-400" />
        <span className="text-sm font-semibold text-gray-200">Controller Designer</span>

        {/* Selector */}
        <div className="relative ml-4">
          <button
            onClick={() => setDropdownOpen(!dropdownOpen)}
            className="flex items-center gap-2 px-3 py-1.5 bg-gray-800 border border-gray-700
                       rounded text-xs hover:bg-gray-750 transition-colors min-w-[200px]"
          >
            <span className="truncate">{controller?.name || 'Select Controller'}</span>
            <ChevronDown size={14} className="ml-auto text-gray-500" />
          </button>
          {dropdownOpen && (
            <>
              <div className="fixed inset-0 z-10" onClick={() => setDropdownOpen(false)} />
              <div className="absolute top-full left-0 mt-1 w-64 bg-gray-800 border border-gray-700
                              rounded-lg shadow-xl z-20 max-h-60 overflow-auto">
                {controllers.map((c) => (
                  <button
                    key={c.id}
                    onClick={() => {
                      setSelectedId(c.id);
                      setDropdownOpen(false);
                    }}
                    className={clsx(
                      'w-full text-left px-3 py-2 text-xs hover:bg-gray-700 transition-colors',
                      c.id === selectedId && 'bg-purple-900/40 text-purple-300'
                    )}
                  >
                    <div className="font-medium">{c.name}</div>
                    <div className="text-gray-500">v{c.version} &middot; {c.controller_type}</div>
                  </button>
                ))}
                {controllers.length === 0 && (
                  <div className="px-3 py-4 text-xs text-gray-500 text-center">No controllers yet</div>
                )}
              </div>
            </>
          )}
        </div>

        <button
          onClick={handleNewController}
          className="flex items-center gap-1 px-3 py-1.5 text-xs bg-gray-800
                     border border-gray-700 rounded hover:bg-gray-750 transition-colors"
        >
          <Plus size={14} />
          New Controller
        </button>

        {controller && (
          <span className="ml-2 px-2 py-0.5 bg-purple-900/40 text-purple-300 text-xs rounded-full font-medium">
            v{controller.version}
          </span>
        )}

        <div className="ml-auto flex items-center gap-2">
          {isDirty && <span className="text-xs text-amber-400 mr-2">Unsaved changes</span>}
          <button
            onClick={handleRevert}
            disabled={!isDirty}
            className="flex items-center gap-1 px-3 py-1.5 text-xs bg-gray-800
                       border border-gray-700 rounded hover:bg-gray-750
                       disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
          >
            <RotateCcw size={14} />
            Revert
          </button>
          <button
            onClick={handleSave}
            disabled={saving || !controller}
            className="flex items-center gap-1 px-4 py-1.5 text-xs bg-purple-600
                       text-white rounded hover:bg-purple-700
                       disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
          >
            {saving ? <Loader2 size={14} className="animate-spin" /> : <Save size={14} />}
            Save
          </button>
        </div>
      </div>

      {/* Error */}
      {error && (
        <div className="flex items-center gap-2 px-4 py-2 bg-red-900/30 border-b border-red-800 text-xs text-red-300">
          <AlertCircle size={14} />
          {error}
          <button onClick={() => setError(null)} className="ml-auto text-red-400 hover:text-red-200">
            Dismiss
          </button>
        </div>
      )}

      {/* Main: 2 panels */}
      {controller ? (
        <div className="flex-1 flex overflow-hidden">
          {/* LEFT PANEL: Parameters (~50%) */}
          <div className="w-1/2 flex flex-col border-r border-gray-800 overflow-hidden">
            <div className="p-3 space-y-2 overflow-auto flex-1">
              {/* Controller name */}
              <div className="flex items-center gap-2">
                <label className="text-xs text-gray-400 w-20 shrink-0">Name</label>
                <input
                  type="text"
                  value={controller.name}
                  onChange={(e) => updateField('name', e.target.value)}
                  className="flex-1 px-2 py-1.5 bg-gray-800 border border-gray-700 rounded
                             text-sm text-gray-200 focus:border-purple-500
                             focus:outline-none focus:ring-1 focus:ring-purple-500/50"
                />
              </div>

              {/* Controller type */}
              <div className="flex items-center gap-2">
                <label className="text-xs text-gray-400 w-20 shrink-0">Type</label>
                <select
                  value={controller.controller_type}
                  onChange={(e) => updateField('controller_type', e.target.value)}
                  className="flex-1 px-2 py-1.5 bg-gray-800 border border-gray-700 rounded
                             text-sm text-gray-200 focus:border-purple-500
                             focus:outline-none"
                >
                  {CONTROLLER_TYPES.map((ct) => (
                    <option key={ct.value} value={ct.value}>
                      {ct.label}
                    </option>
                  ))}
                </select>
              </div>

              {/* Pitch Control Section */}
              <SectionHeader title="Pitch Control" icon={Settings} />
              <div className="bg-gray-900/50 rounded-lg p-3 space-y-2 border border-gray-800">
                <ParamInput
                  label="PCMode"
                  value={controller.pcmode}
                  onChange={(v) => updateField('pcmode', v)}
                />
                <ParamInput
                  label="TPCOn"
                  value={params.tpcon ?? 0}
                  unit="s"
                  onChange={(v) => updateParam('tpcon', v)}
                />
                <ParamInput
                  label="Pitch Min"
                  value={params.pitch_min ?? 0}
                  unit="deg"
                  onChange={(v) => updateParam('pitch_min', v)}
                />
                <ParamInput
                  label="Pitch Max"
                  value={params.pitch_max ?? 90}
                  unit="deg"
                  onChange={(v) => updateParam('pitch_max', v)}
                />
                <ParamInput
                  label="Pitch Rate Limit"
                  value={params.pitch_rate_limit ?? 8}
                  unit="deg/s"
                  onChange={(v) => updateParam('pitch_rate_limit', v)}
                />
                <ParamInput
                  label="KP (Pitch Gain)"
                  value={params.kp_pitch ?? 0}
                  unit="s"
                  onChange={(v) => updateParam('kp_pitch', v)}
                />
                <ParamInput
                  label="KI (Pitch Gain)"
                  value={params.ki_pitch ?? 0}
                  unit="-"
                  onChange={(v) => updateParam('ki_pitch', v)}
                />
              </div>

              {/* Torque Control Section */}
              <SectionHeader title="Torque Control" icon={Settings} />
              <div className="bg-gray-900/50 rounded-lg p-3 space-y-2 border border-gray-800">
                <ParamInput
                  label="VSContrl"
                  value={controller.vscontrl}
                  onChange={(v) => updateField('vscontrl', v)}
                />
                <ParamInput
                  label="Rated Gen Speed"
                  value={params.rated_gen_speed ?? 0}
                  unit="rpm"
                  onChange={(v) => updateParam('rated_gen_speed', v)}
                />
                <ParamInput
                  label="Rgn2K"
                  value={params.rgn2k ?? 0}
                  unit="N-m/rpm^2"
                  onChange={(v) => updateParam('rgn2k', v)}
                />
                <ParamInput
                  label="Rated Torque"
                  value={params.rated_torque ?? 0}
                  unit="N-m"
                  onChange={(v) => updateParam('rated_torque', v)}
                />
                <ParamInput
                  label="Min Gen Torque"
                  value={params.min_gen_torque ?? 0}
                  unit="N-m"
                  onChange={(v) => updateParam('min_gen_torque', v)}
                />
                <ParamInput
                  label="Max Gen Torque"
                  value={params.max_gen_torque ?? 0}
                  unit="N-m"
                  onChange={(v) => updateParam('max_gen_torque', v)}
                />
              </div>

              {/* Other Section */}
              <SectionHeader title="Other Parameters" icon={Settings} />
              <div className="bg-gray-900/50 rounded-lg p-3 space-y-2 border border-gray-800">
                <ParamInput
                  label="GenModel"
                  value={params.gen_model ?? 1}
                  onChange={(v) => updateParam('gen_model', v)}
                />
                <ParamInput
                  label="Generator Efficiency"
                  value={params.gen_eff ?? 94.4}
                  unit="%"
                  onChange={(v) => updateParam('gen_eff', v)}
                />
                <ParamInput
                  label="DLL Filename"
                  value={controller.dll_filename ?? ''}
                  type="text"
                  onChange={(v) => updateField('dll_filename', v || null)}
                />
                <ParamInput
                  label="DLL Procedure"
                  value={controller.dll_procname ?? ''}
                  type="text"
                  onChange={(v) => updateField('dll_procname', v || null)}
                />
              </div>

              {/* ROSCO Gain Schedule */}
              {params.gain_schedule && Array.isArray(params.gain_schedule) && (
                <>
                  <SectionHeader title="Gain Schedule" icon={Settings} />
                  <GainScheduleTable
                    schedule={params.gain_schedule}
                    onChange={(s) => updateParam('gain_schedule', s)}
                  />
                </>
              )}
            </div>
          </div>

          {/* RIGHT PANEL: Visualization (~50%) */}
          <div className="w-1/2 flex flex-col overflow-hidden">
            {/* Tabs */}
            <div className="flex border-b border-gray-800 shrink-0">
              {(['regions', 'pitch_schedule', 'preview'] as const).map((tab) => (
                <button
                  key={tab}
                  onClick={() => setRightTab(tab)}
                  className={clsx(
                    'flex-1 px-3 py-2 text-xs font-medium transition-colors',
                    rightTab === tab
                      ? 'bg-gray-800 text-purple-400 border-b-2 border-purple-400'
                      : 'text-gray-400 hover:text-gray-200 hover:bg-gray-900'
                  )}
                >
                  {tab === 'regions'
                    ? 'Operating Regions'
                    : tab === 'pitch_schedule'
                    ? 'Pitch Schedule'
                    : 'File Preview'}
                </button>
              ))}
            </div>

            <div className="flex-1 overflow-auto">
              {rightTab === 'regions' ? (
                <div className="p-4">
                  <TorqueSpeedPlot controller={controller} />
                  <div className="mt-4 grid grid-cols-2 gap-3">
                    <div className="bg-gray-900/50 rounded-lg p-3 border border-gray-800">
                      <div className="text-xs text-gray-500 mb-1">Rated Power</div>
                      <div className="text-sm font-mono text-gray-200">
                        {(
                          ((params.rated_torque ?? 0) * (params.rated_gen_speed ?? 0) * Math.PI) /
                          30 /
                          1e6
                        ).toFixed(2)}{' '}
                        MW
                      </div>
                    </div>
                    <div className="bg-gray-900/50 rounded-lg p-3 border border-gray-800">
                      <div className="text-xs text-gray-500 mb-1">Rated Speed</div>
                      <div className="text-sm font-mono text-gray-200">
                        {(params.rated_gen_speed ?? 0).toFixed(1)} rpm
                      </div>
                    </div>
                    <div className="bg-gray-900/50 rounded-lg p-3 border border-gray-800">
                      <div className="text-xs text-gray-500 mb-1">Rated Torque</div>
                      <div className="text-sm font-mono text-gray-200">
                        {(params.rated_torque ?? 0).toFixed(1)} N-m
                      </div>
                    </div>
                    <div className="bg-gray-900/50 rounded-lg p-3 border border-gray-800">
                      <div className="text-xs text-gray-500 mb-1">Rgn2K</div>
                      <div className="text-sm font-mono text-gray-200">
                        {(params.rgn2k ?? 0).toFixed(7)} N-m/rpm^2
                      </div>
                    </div>
                  </div>
                </div>
              ) : rightTab === 'pitch_schedule' ? (
                <div className="p-4">
                  <PlotPanel
                    data={gainSchedulePlot.data}
                    layout={gainSchedulePlot.layout}
                    title="Gain Schedule (KP & KI vs Wind Speed)"
                    height={350}
                  />
                </div>
              ) : (
                <div className="p-2 h-full">
                  <FilePreview
                    content={filePreview}
                    title="ServoDyn Input File"
                    isLoading={previewLoading}
                    className="h-full"
                  />
                </div>
              )}
            </div>
          </div>
        </div>
      ) : (
        <div className="flex-1 flex items-center justify-center">
          <div className="text-center text-gray-500">
            <Cpu size={64} className="mx-auto mb-4 opacity-20" />
            <p className="text-lg font-medium">No controller selected</p>
            <p className="text-sm mt-2">Select a controller or create a new one</p>
            <button
              onClick={handleNewController}
              className="mt-4 px-4 py-2 bg-purple-600 text-white rounded
                         hover:bg-purple-700 transition-colors text-sm"
            >
              Create New Controller
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
