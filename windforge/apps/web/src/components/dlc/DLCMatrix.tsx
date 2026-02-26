import { useState, useEffect, useMemo, useCallback } from 'react';
import { useParams } from 'react-router-dom';
import {
  Save,
  FolderOpen,
  ChevronDown,
  ChevronRight,
  Zap,
  Shield,
  Wind,
  Settings2,
  Hash,
  RotateCcw,
} from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';
import apiClient from '@/api/client';
import StatusBadge from '@/components/common/StatusBadge';

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface DLCCaseSpec {
  dlc_number: string;
  wind_speeds: number[];
  seeds: number;
  yaw_misalignments: number[];
}

interface TurbSimParams {
  turbulence_model: string;
  iec_standard: string;
  iec_turbc: string;
  grid_height: number;
  grid_width: number;
  num_grid_z: number;
  num_grid_y: number;
  time_step: number;
  analysis_time: number;
  ref_height: number;
}

interface DLCDefinition {
  id: string;
  project_id: string;
  turbine_model_id: string;
  name: string;
  dlc_cases: DLCCaseSpec[] | null;
  turbsim_params: TurbSimParams | null;
  total_case_count: number;
  status: string;
  created_at: string;
}

interface TurbineModel {
  id: string;
  name: string;
}

// ---------------------------------------------------------------------------
// DLC metadata
// ---------------------------------------------------------------------------

interface DLCMeta {
  number: string;
  description: string;
  windModel: string;
  analysisType: 'fatigue' | 'ultimate';
}

const DLC_CATALOG: DLCMeta[] = [
  { number: '1.1', description: 'Power production', windModel: 'NTM', analysisType: 'fatigue' },
  { number: '1.2', description: 'Power production (no turb. SF)', windModel: 'NTM', analysisType: 'fatigue' },
  { number: '1.3', description: 'Power production', windModel: 'ETM', analysisType: 'ultimate' },
  { number: '1.4', description: 'Power production', windModel: 'ECD', analysisType: 'ultimate' },
  { number: '1.5', description: 'Power production', windModel: 'EWS', analysisType: 'ultimate' },
  { number: '2.1', description: 'Power prod + fault', windModel: 'NTM', analysisType: 'ultimate' },
  { number: '2.2', description: 'Power prod + fault', windModel: 'NTM', analysisType: 'ultimate' },
  { number: '2.3', description: 'Power prod + EOG', windModel: 'EOG', analysisType: 'ultimate' },
  { number: '3.1', description: 'Start-up', windModel: 'NWP+EOG', analysisType: 'fatigue' },
  { number: '3.2', description: 'Start-up', windModel: 'EDC', analysisType: 'ultimate' },
  { number: '3.3', description: 'Start-up', windModel: 'ECD', analysisType: 'ultimate' },
  { number: '4.1', description: 'Normal shutdown', windModel: 'NWP+EOG', analysisType: 'fatigue' },
  { number: '4.2', description: 'Normal shutdown', windModel: 'NWP+EOG', analysisType: 'ultimate' },
  { number: '5.1', description: 'Emergency shutdown', windModel: 'NTM', analysisType: 'ultimate' },
  { number: '6.1', description: 'Parked (standstill)', windModel: 'EWM 50-yr', analysisType: 'ultimate' },
  { number: '6.2', description: 'Parked (standstill)', windModel: 'EWM 50-yr', analysisType: 'ultimate' },
  { number: '6.3', description: 'Parked (standstill)', windModel: 'EWM 1-yr', analysisType: 'ultimate' },
  { number: '6.4', description: 'Parked (standstill)', windModel: 'NTM', analysisType: 'fatigue' },
  { number: '7.1', description: 'Parked + fault', windModel: 'EWM 1-yr', analysisType: 'ultimate' },
  { number: '8.1', description: 'Transport/installation', windModel: 'NTM', analysisType: 'ultimate' },
];

const DEFAULT_YAW_OPTIONS = [-8, -4, 0, 4, 8];

// ---------------------------------------------------------------------------
// Per-DLC row state
// ---------------------------------------------------------------------------

interface DLCRowState {
  enabled: boolean;
  expanded: boolean;
  windSpeedMin: number;
  windSpeedMax: number;
  windSpeedStep: number;
  seeds: number;
  yawMisalignments: number[];
}

function defaultRowState(): DLCRowState {
  return {
    enabled: false,
    expanded: false,
    windSpeedMin: 4,
    windSpeedMax: 24,
    windSpeedStep: 2,
    seeds: 6,
    yawMisalignments: [0],
  };
}

function generateWindSpeeds(min: number, max: number, step: number): number[] {
  const speeds: number[] = [];
  for (let v = min; v <= max + 0.001; v += step) {
    speeds.push(parseFloat(v.toFixed(1)));
  }
  return speeds;
}

function countCases(row: DLCRowState): number {
  if (!row.enabled) return 0;
  const speeds = generateWindSpeeds(row.windSpeedMin, row.windSpeedMax, row.windSpeedStep);
  return speeds.length * row.seeds * row.yawMisalignments.length;
}

// ---------------------------------------------------------------------------
// Default TurbSim params
// ---------------------------------------------------------------------------

function defaultTurbSimParams(): TurbSimParams {
  return {
    turbulence_model: 'IECKAI',
    iec_standard: '1-ED3',
    iec_turbc: 'B',
    grid_height: 180,
    grid_width: 180,
    num_grid_z: 25,
    num_grid_y: 25,
    time_step: 0.05,
    analysis_time: 660,
    ref_height: 119,
  };
}

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------

export default function DLCMatrix() {
  const { projectId } = useParams<{ projectId: string }>();

  // ---- State ----
  const [name, setName] = useState('Default DLC Matrix');
  const [turbineModels, setTurbineModels] = useState<TurbineModel[]>([]);
  const [selectedTurbineId, setSelectedTurbineId] = useState('');
  const [existingDefinitions, setExistingDefinitions] = useState<DLCDefinition[]>([]);
  const [activeDefinitionId, setActiveDefinitionId] = useState<string | null>(null);
  const [saving, setSaving] = useState(false);

  const [rows, setRows] = useState<Record<string, DLCRowState>>(() => {
    const init: Record<string, DLCRowState> = {};
    DLC_CATALOG.forEach((dlc) => {
      init[dlc.number] = defaultRowState();
    });
    return init;
  });

  const [turbSimParams, setTurbSimParams] = useState<TurbSimParams>(defaultTurbSimParams());

  // ---- Fetch turbine models & existing definitions ----
  useEffect(() => {
    if (!projectId) return;

    apiClient
      .get(`/projects/${projectId}/turbine-models`)
      .then((res) => {
        const models = res.data as TurbineModel[];
        setTurbineModels(models);
        if (models.length > 0 && !selectedTurbineId) {
          setSelectedTurbineId(models[0].id);
        }
      })
      .catch(() => {
        // Turbine models might not exist yet
      });

    apiClient
      .get(`/projects/${projectId}/dlc-definitions`)
      .then((res) => {
        setExistingDefinitions(res.data as DLCDefinition[]);
      })
      .catch(() => {});
  }, [projectId]);

  // ---- Load an existing definition ----
  const loadDefinition = useCallback(
    (def: DLCDefinition) => {
      setName(def.name);
      setSelectedTurbineId(def.turbine_model_id);
      setActiveDefinitionId(def.id);

      if (def.turbsim_params) {
        setTurbSimParams(def.turbsim_params);
      }

      const newRows: Record<string, DLCRowState> = {};
      DLC_CATALOG.forEach((dlc) => {
        newRows[dlc.number] = defaultRowState();
      });

      if (def.dlc_cases) {
        def.dlc_cases.forEach((c) => {
          if (!newRows[c.dlc_number]) return;
          const speeds = [...c.wind_speeds].sort((a, b) => a - b);
          const step =
            speeds.length > 1 ? parseFloat((speeds[1] - speeds[0]).toFixed(1)) : 2;
          newRows[c.dlc_number] = {
            enabled: true,
            expanded: false,
            windSpeedMin: speeds[0] ?? 4,
            windSpeedMax: speeds[speeds.length - 1] ?? 24,
            windSpeedStep: step,
            seeds: c.seeds,
            yawMisalignments: c.yaw_misalignments,
          };
        });
      }

      setRows(newRows);
      toast.success(`Loaded "${def.name}"`);
    },
    [],
  );

  // ---- Row helpers ----
  const updateRow = useCallback((dlcNum: string, patch: Partial<DLCRowState>) => {
    setRows((prev) => ({
      ...prev,
      [dlcNum]: { ...prev[dlcNum], ...patch },
    }));
  }, []);

  const toggleYaw = useCallback((dlcNum: string, yaw: number) => {
    setRows((prev) => {
      const cur = prev[dlcNum];
      const yaws = cur.yawMisalignments.includes(yaw)
        ? cur.yawMisalignments.filter((y) => y !== yaw)
        : [...cur.yawMisalignments, yaw].sort((a, b) => a - b);
      return { ...prev, [dlcNum]: { ...cur, yawMisalignments: yaws.length > 0 ? yaws : [0] } };
    });
  }, []);

  // ---- Derived totals ----
  const caseCountByDLC = useMemo(() => {
    const counts: Record<string, number> = {};
    DLC_CATALOG.forEach((dlc) => {
      counts[dlc.number] = countCases(rows[dlc.number]);
    });
    return counts;
  }, [rows]);

  const totalCases = useMemo(
    () => Object.values(caseCountByDLC).reduce((s, c) => s + c, 0),
    [caseCountByDLC],
  );

  const enabledCount = useMemo(
    () => DLC_CATALOG.filter((d) => rows[d.number].enabled).length,
    [rows],
  );

  // ---- Build API payload ----
  const buildPayload = useCallback(() => {
    const dlcCases: DLCCaseSpec[] = DLC_CATALOG.filter((d) => rows[d.number].enabled).map(
      (d) => {
        const r = rows[d.number];
        return {
          dlc_number: d.number,
          wind_speeds: generateWindSpeeds(r.windSpeedMin, r.windSpeedMax, r.windSpeedStep),
          seeds: r.seeds,
          yaw_misalignments: r.yawMisalignments,
        };
      },
    );
    return {
      name,
      turbine_model_id: selectedTurbineId,
      dlc_cases: dlcCases,
      turbsim_params: turbSimParams,
    };
  }, [name, selectedTurbineId, rows, turbSimParams]);

  // ---- Save ----
  const handleSave = useCallback(async () => {
    if (!projectId) return;
    setSaving(true);
    try {
      const payload = buildPayload();
      if (activeDefinitionId) {
        await apiClient.put(
          `/projects/${projectId}/dlc-definitions/${activeDefinitionId}`,
          payload,
        );
        toast.success('DLC definition updated');
      } else {
        const res = await apiClient.post(
          `/projects/${projectId}/dlc-definitions`,
          payload,
        );
        setActiveDefinitionId(res.data.id);
        toast.success('DLC definition created');
      }
      // Refresh list
      const res = await apiClient.get(`/projects/${projectId}/dlc-definitions`);
      setExistingDefinitions(res.data as DLCDefinition[]);
    } catch (err: any) {
      toast.error(err?.response?.data?.detail ?? 'Failed to save');
    } finally {
      setSaving(false);
    }
  }, [projectId, activeDefinitionId, buildPayload]);

  // ---- Generate cases (save + trigger) ----
  const handleGenerateCases = useCallback(async () => {
    if (enabledCount === 0) {
      toast.error('Enable at least one DLC');
      return;
    }
    await handleSave();
    toast.success(`${totalCases} cases ready for simulation`);
  }, [enabledCount, totalCases, handleSave]);

  // ---- TurbSim param helpers ----
  const updateTurbSim = useCallback((patch: Partial<TurbSimParams>) => {
    setTurbSimParams((prev) => ({ ...prev, ...patch }));
  }, []);

  // ---- Render ----
  return (
    <div className="flex h-full gap-6 overflow-hidden">
      {/* ====== Main area ====== */}
      <div className="flex flex-1 flex-col overflow-hidden">
        {/* Top bar */}
        <div className="mb-4 flex flex-wrap items-center gap-3">
          <input
            type="text"
            value={name}
            onChange={(e) => setName(e.target.value)}
            className="input-field max-w-xs"
            placeholder="DLC Definition Name"
          />

          <select
            value={selectedTurbineId}
            onChange={(e) => setSelectedTurbineId(e.target.value)}
            className="input-field max-w-[200px]"
          >
            <option value="">Select Turbine</option>
            {turbineModels.map((m) => (
              <option key={m.id} value={m.id}>
                {m.name}
              </option>
            ))}
          </select>

          <span className="inline-flex items-center gap-1.5 rounded-full bg-accent-500/20 px-3 py-1 text-xs font-semibold text-accent-300 ring-1 ring-inset ring-accent-500/30">
            <Hash size={12} />
            {totalCases.toLocaleString()} cases
          </span>

          <div className="ml-auto flex gap-2">
            {/* Load dropdown */}
            {existingDefinitions.length > 0 && (
              <div className="relative group">
                <button className="btn-secondary text-xs">
                  <FolderOpen size={14} />
                  Load
                </button>
                <div className="invisible absolute right-0 z-20 mt-1 w-64 rounded-lg border border-slate-600 bg-surface-dark-secondary shadow-xl group-hover:visible">
                  {existingDefinitions.map((def) => (
                    <button
                      key={def.id}
                      onClick={() => loadDefinition(def)}
                      className="flex w-full items-center justify-between px-4 py-2 text-left text-sm text-slate-200 hover:bg-surface-dark-tertiary first:rounded-t-lg last:rounded-b-lg"
                    >
                      <span>{def.name}</span>
                      <StatusBadge status={def.status} size="sm" />
                    </button>
                  ))}
                </div>
              </div>
            )}

            <button
              onClick={handleSave}
              disabled={saving}
              className="btn-secondary text-xs"
            >
              <Save size={14} />
              {saving ? 'Saving...' : 'Save'}
            </button>
          </div>
        </div>

        {/* DLC Grid */}
        <div className="flex-1 overflow-auto rounded-xl border border-slate-700/50 bg-surface-dark-secondary">
          <table className="w-full text-sm">
            <thead className="sticky top-0 z-10 bg-surface-dark-secondary">
              <tr className="border-b border-slate-700 text-slate-400">
                <th className="w-8 px-3 py-3" />
                <th className="w-10 px-3 py-3 text-left font-semibold">En</th>
                <th className="px-3 py-3 text-left font-semibold">DLC</th>
                <th className="px-3 py-3 text-left font-semibold">Description</th>
                <th className="px-3 py-3 text-left font-semibold">Wind Model</th>
                <th className="px-3 py-3 text-left font-semibold">Type</th>
                <th className="px-3 py-3 text-right font-semibold">Cases</th>
              </tr>
            </thead>
            <tbody>
              {DLC_CATALOG.map((dlc) => {
                const row = rows[dlc.number];
                const cases = caseCountByDLC[dlc.number];
                const isFatigue = dlc.analysisType === 'fatigue';

                return (
                  <DLCRow
                    key={dlc.number}
                    dlc={dlc}
                    row={row}
                    cases={cases}
                    isFatigue={isFatigue}
                    onToggleEnabled={() => updateRow(dlc.number, { enabled: !row.enabled })}
                    onToggleExpand={() => updateRow(dlc.number, { expanded: !row.expanded })}
                    onUpdateRow={(patch) => updateRow(dlc.number, patch)}
                    onToggleYaw={(yaw) => toggleYaw(dlc.number, yaw)}
                  />
                );
              })}
            </tbody>
          </table>
        </div>
      </div>

      {/* ====== Right sidebar ====== */}
      <div className="w-80 shrink-0 overflow-y-auto rounded-xl border border-slate-700/50 bg-surface-dark-secondary p-5">
        {/* TurbSim Parameters */}
        <div className="mb-6">
          <h3 className="mb-4 flex items-center gap-2 text-sm font-semibold text-slate-100">
            <Settings2 size={16} className="text-accent-400" />
            TurbSim Parameters
          </h3>

          <div className="space-y-3">
            <div>
              <label className="label">Turbulence Model</label>
              <select
                value={turbSimParams.turbulence_model}
                onChange={(e) => updateTurbSim({ turbulence_model: e.target.value })}
                className="input-field"
              >
                <option value="IECKAI">IECKAI (Kaimal)</option>
                <option value="IECVKM">IECVKM (von Karman)</option>
                <option value="GP_LLJ">GP_LLJ (Great Plains LLJ)</option>
              </select>
            </div>

            <div className="grid grid-cols-2 gap-3">
              <div>
                <label className="label">IEC Standard</label>
                <select
                  value={turbSimParams.iec_standard}
                  onChange={(e) => updateTurbSim({ iec_standard: e.target.value })}
                  className="input-field"
                >
                  <option value="1-ED3">IEC 61400-1 Ed.3</option>
                  <option value="1-ED4">IEC 61400-1 Ed.4</option>
                </select>
              </div>
              <div>
                <label className="label">Turb. Category</label>
                <select
                  value={turbSimParams.iec_turbc}
                  onChange={(e) => updateTurbSim({ iec_turbc: e.target.value })}
                  className="input-field"
                >
                  <option value="A">A (high)</option>
                  <option value="B">B (medium)</option>
                  <option value="C">C (low)</option>
                </select>
              </div>
            </div>

            <div className="grid grid-cols-2 gap-3">
              <TurbSimField
                label="Grid Height"
                unit="m"
                value={turbSimParams.grid_height}
                onChange={(v) => updateTurbSim({ grid_height: v })}
              />
              <TurbSimField
                label="Grid Width"
                unit="m"
                value={turbSimParams.grid_width}
                onChange={(v) => updateTurbSim({ grid_width: v })}
              />
            </div>

            <div className="grid grid-cols-2 gap-3">
              <TurbSimField
                label="Z Points"
                value={turbSimParams.num_grid_z}
                onChange={(v) => updateTurbSim({ num_grid_z: v })}
                step={1}
                min={3}
              />
              <TurbSimField
                label="Y Points"
                value={turbSimParams.num_grid_y}
                onChange={(v) => updateTurbSim({ num_grid_y: v })}
                step={1}
                min={3}
              />
            </div>

            <div className="grid grid-cols-2 gap-3">
              <TurbSimField
                label="Time Step"
                unit="s"
                value={turbSimParams.time_step}
                onChange={(v) => updateTurbSim({ time_step: v })}
                step={0.01}
                min={0.01}
              />
              <TurbSimField
                label="Analysis Time"
                unit="s"
                value={turbSimParams.analysis_time}
                onChange={(v) => updateTurbSim({ analysis_time: v })}
                step={10}
                min={60}
              />
            </div>

            <TurbSimField
              label="Reference Height"
              unit="m"
              value={turbSimParams.ref_height}
              onChange={(v) => updateTurbSim({ ref_height: v })}
              min={10}
            />
          </div>
        </div>

        {/* Divider */}
        <div className="my-5 border-t border-slate-700" />

        {/* Case count breakdown */}
        <div className="mb-6">
          <h3 className="mb-3 text-sm font-semibold text-slate-100">Case Breakdown</h3>
          <div className="space-y-1.5 max-h-64 overflow-y-auto">
            {DLC_CATALOG.filter((d) => rows[d.number].enabled).map((dlc) => (
              <div
                key={dlc.number}
                className="flex items-center justify-between rounded-md bg-surface-dark px-3 py-1.5 text-xs"
              >
                <span className="font-mono text-slate-300">DLC {dlc.number}</span>
                <span className="font-semibold text-slate-100">
                  {caseCountByDLC[dlc.number].toLocaleString()}
                </span>
              </div>
            ))}
            {enabledCount === 0 && (
              <p className="py-4 text-center text-xs text-slate-500">
                No DLCs enabled
              </p>
            )}
          </div>
        </div>

        {/* Divider */}
        <div className="my-5 border-t border-slate-700" />

        {/* Summary + action */}
        <div className="space-y-3">
          <div className="flex items-center justify-between text-sm">
            <span className="text-slate-400">Enabled DLCs</span>
            <span className="font-semibold text-slate-100">{enabledCount}</span>
          </div>
          <div className="flex items-center justify-between text-sm">
            <span className="text-slate-400">Total Cases</span>
            <span className="font-semibold text-accent-300">{totalCases.toLocaleString()}</span>
          </div>
          <button
            onClick={handleGenerateCases}
            disabled={enabledCount === 0 || saving}
            className="btn-primary w-full"
          >
            <Zap size={16} />
            Generate Cases
          </button>
          <button
            onClick={() => {
              const init: Record<string, DLCRowState> = {};
              DLC_CATALOG.forEach((d) => {
                init[d.number] = defaultRowState();
              });
              setRows(init);
              setTurbSimParams(defaultTurbSimParams());
              setActiveDefinitionId(null);
              toast('Reset to defaults');
            }}
            className="btn-secondary w-full text-xs"
          >
            <RotateCcw size={14} />
            Reset All
          </button>
        </div>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Sub-components
// ---------------------------------------------------------------------------

function TurbSimField({
  label,
  unit,
  value,
  onChange,
  step = 1,
  min = 0,
}: {
  label: string;
  unit?: string;
  value: number;
  onChange: (v: number) => void;
  step?: number;
  min?: number;
}) {
  return (
    <div>
      <label className="label">
        {label}
        {unit && <span className="ml-1 text-slate-500">[{unit}]</span>}
      </label>
      <input
        type="number"
        value={value}
        onChange={(e) => {
          const v = parseFloat(e.target.value);
          if (!isNaN(v) && v >= min) onChange(v);
        }}
        step={step}
        min={min}
        className="input-field font-mono text-xs"
      />
    </div>
  );
}

interface DLCRowProps {
  dlc: DLCMeta;
  row: DLCRowState;
  cases: number;
  isFatigue: boolean;
  onToggleEnabled: () => void;
  onToggleExpand: () => void;
  onUpdateRow: (patch: Partial<DLCRowState>) => void;
  onToggleYaw: (yaw: number) => void;
}

function DLCRow({
  dlc,
  row,
  cases,
  isFatigue,
  onToggleEnabled,
  onToggleExpand,
  onUpdateRow,
  onToggleYaw,
}: DLCRowProps) {
  return (
    <>
      {/* Main row */}
      <tr
        className={clsx(
          'border-b border-slate-700/50 transition-colors',
          row.enabled
            ? 'bg-surface-dark hover:bg-surface-dark-tertiary/30'
            : 'bg-surface-dark-secondary/50 opacity-60 hover:opacity-80',
          row.enabled && isFatigue && 'border-l-2 border-l-emerald-500/50',
          row.enabled && !isFatigue && 'border-l-2 border-l-red-500/50',
        )}
      >
        {/* Expand */}
        <td className="px-3 py-2.5">
          {row.enabled && (
            <button
              onClick={onToggleExpand}
              className="text-slate-400 hover:text-slate-200 transition-colors"
            >
              {row.expanded ? <ChevronDown size={14} /> : <ChevronRight size={14} />}
            </button>
          )}
        </td>

        {/* Enable */}
        <td className="px-3 py-2.5">
          <input
            type="checkbox"
            checked={row.enabled}
            onChange={onToggleEnabled}
            className="h-4 w-4 rounded border-slate-500 bg-surface-dark text-accent-500 focus:ring-accent-500 focus:ring-offset-0"
          />
        </td>

        {/* DLC number */}
        <td className="px-3 py-2.5 font-mono font-semibold text-slate-100">
          {dlc.number}
        </td>

        {/* Description */}
        <td className="px-3 py-2.5 text-slate-300">{dlc.description}</td>

        {/* Wind model badge */}
        <td className="px-3 py-2.5">
          <span className="inline-flex items-center gap-1 rounded-md bg-blue-500/15 px-2 py-0.5 text-xs font-medium text-blue-300 ring-1 ring-inset ring-blue-500/25">
            <Wind size={11} />
            {dlc.windModel}
          </span>
        </td>

        {/* Analysis type badge */}
        <td className="px-3 py-2.5">
          <span
            className={clsx(
              'inline-flex items-center gap-1 rounded-md px-2 py-0.5 text-xs font-medium ring-1 ring-inset',
              isFatigue
                ? 'bg-emerald-500/15 text-emerald-300 ring-emerald-500/25'
                : 'bg-red-500/15 text-red-300 ring-red-500/25',
            )}
          >
            {isFatigue ? <RotateCcw size={11} /> : <Shield size={11} />}
            {dlc.analysisType}
          </span>
        </td>

        {/* Case count */}
        <td className="px-3 py-2.5 text-right font-mono text-slate-300">
          {row.enabled ? cases.toLocaleString() : '-'}
        </td>
      </tr>

      {/* Expanded config */}
      {row.enabled && row.expanded && (
        <tr className="border-b border-slate-700/50 bg-surface-dark/80">
          <td colSpan={7} className="px-6 py-4">
            <div className="flex flex-wrap gap-6">
              {/* Wind speed range */}
              <div>
                <p className="mb-2 text-xs font-semibold text-slate-400 uppercase tracking-wider">
                  Wind Speed Range
                </p>
                <div className="flex items-center gap-2">
                  <div>
                    <label className="text-[10px] text-slate-500">Vin (m/s)</label>
                    <input
                      type="number"
                      value={row.windSpeedMin}
                      onChange={(e) => {
                        const v = parseFloat(e.target.value);
                        if (!isNaN(v) && v >= 0) onUpdateRow({ windSpeedMin: v });
                      }}
                      step={1}
                      min={0}
                      className="input-field w-20 text-xs font-mono"
                    />
                  </div>
                  <span className="mt-4 text-slate-500">to</span>
                  <div>
                    <label className="text-[10px] text-slate-500">Vout (m/s)</label>
                    <input
                      type="number"
                      value={row.windSpeedMax}
                      onChange={(e) => {
                        const v = parseFloat(e.target.value);
                        if (!isNaN(v) && v > row.windSpeedMin) onUpdateRow({ windSpeedMax: v });
                      }}
                      step={1}
                      min={row.windSpeedMin + 1}
                      className="input-field w-20 text-xs font-mono"
                    />
                  </div>
                  <div>
                    <label className="text-[10px] text-slate-500">Step (m/s)</label>
                    <input
                      type="number"
                      value={row.windSpeedStep}
                      onChange={(e) => {
                        const v = parseFloat(e.target.value);
                        if (!isNaN(v) && v > 0) onUpdateRow({ windSpeedStep: v });
                      }}
                      step={0.5}
                      min={0.5}
                      className="input-field w-20 text-xs font-mono"
                    />
                  </div>
                </div>
                <p className="mt-1 text-[10px] text-slate-500">
                  {generateWindSpeeds(row.windSpeedMin, row.windSpeedMax, row.windSpeedStep).length}{' '}
                  wind speeds:{' '}
                  {generateWindSpeeds(row.windSpeedMin, row.windSpeedMax, row.windSpeedStep)
                    .map((v) => v.toFixed(0))
                    .join(', ')}
                </p>
              </div>

              {/* Seeds */}
              <div>
                <p className="mb-2 text-xs font-semibold text-slate-400 uppercase tracking-wider">
                  Seeds
                </p>
                <div className="flex items-center gap-2">
                  <button
                    onClick={() => onUpdateRow({ seeds: Math.max(1, row.seeds - 1) })}
                    className="flex h-8 w-8 items-center justify-center rounded-md border border-slate-600 bg-surface-dark-tertiary text-slate-300 hover:bg-slate-600 transition-colors"
                  >
                    -
                  </button>
                  <span className="w-8 text-center font-mono font-semibold text-slate-100">
                    {row.seeds}
                  </span>
                  <button
                    onClick={() => onUpdateRow({ seeds: Math.min(30, row.seeds + 1) })}
                    className="flex h-8 w-8 items-center justify-center rounded-md border border-slate-600 bg-surface-dark-tertiary text-slate-300 hover:bg-slate-600 transition-colors"
                  >
                    +
                  </button>
                </div>
              </div>

              {/* Yaw misalignment */}
              <div>
                <p className="mb-2 text-xs font-semibold text-slate-400 uppercase tracking-wider">
                  Yaw Misalignment
                </p>
                <div className="flex flex-wrap gap-1.5">
                  {DEFAULT_YAW_OPTIONS.map((yaw) => (
                    <button
                      key={yaw}
                      onClick={() => onToggleYaw(yaw)}
                      className={clsx(
                        'rounded-md px-2.5 py-1 text-xs font-mono font-medium transition-colors ring-1 ring-inset',
                        row.yawMisalignments.includes(yaw)
                          ? 'bg-accent-500/20 text-accent-300 ring-accent-500/40'
                          : 'bg-surface-dark-tertiary text-slate-400 ring-slate-600 hover:ring-slate-500',
                      )}
                    >
                      {yaw > 0 ? `+${yaw}` : yaw}°
                    </button>
                  ))}
                </div>
              </div>
            </div>
          </td>
        </tr>
      )}
    </>
  );
}
