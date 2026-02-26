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
  CheckCircle2,
  XCircle,
  Zap,
  ShieldCheck,
} from 'lucide-react';
import clsx from 'clsx';

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface Tower {
  id: string;
  name: string;
  version: number;
  tower_height: number;
}

interface Blade {
  id: string;
  name: string;
  version: number;
  blade_length: number;
}

interface Controller {
  id: string;
  name: string;
  version: number;
  controller_type: string;
}

interface TurbineModel {
  id: string;
  project_id: string;
  name: string;
  version: number;
  tower_id: string | null;
  blade_id: string | null;
  controller_id: string | null;
  gearbox_ratio: number | null;
  generator_inertia: number | null;
  drivetrain_stiffness: number | null;
  drivetrain_damping: number | null;
  hub_mass: number | null;
  hub_inertia: number | null;
  nacelle_mass: number | null;
  nacelle_inertia: number | null;
  overhang: number | null;
  shaft_tilt: number | null;
  precone: number | null;
  rotor_speed_rated: number | null;
  dof_flags: Record<string, boolean> | null;
  is_active: boolean;
  created_at: string;
}

const API_BASE = '/api/v1';

const DOF_FLAGS = [
  { key: 'FlapDOF1', label: 'Flap DOF 1' },
  { key: 'FlapDOF2', label: 'Flap DOF 2' },
  { key: 'EdgeDOF1', label: 'Edge DOF 1' },
  { key: 'TwrFADOF1', label: 'Tower FA DOF 1' },
  { key: 'TwrFADOF2', label: 'Tower FA DOF 2' },
  { key: 'TwrSSDOF1', label: 'Tower SS DOF 1' },
  { key: 'TwrSSDOF2', label: 'Tower SS DOF 2' },
  { key: 'DrTrDOF', label: 'Drivetrain DOF' },
  { key: 'GenDOF', label: 'Generator DOF' },
  { key: 'YawDOF', label: 'Yaw DOF' },
  { key: 'PtfmSgDOF', label: 'Platform Surge DOF' },
  { key: 'PtfmSwDOF', label: 'Platform Sway DOF' },
  { key: 'PtfmHvDOF', label: 'Platform Heave DOF' },
];

// ---------------------------------------------------------------------------
// Helper components
// ---------------------------------------------------------------------------

function ParamInput({
  label,
  value,
  unit,
  onChange,
}: {
  label: string;
  value: number | null;
  unit?: string;
  onChange: (v: number | null) => void;
}) {
  return (
    <div className="flex items-center gap-2">
      <label className="text-xs text-gray-400 w-40 shrink-0 truncate" title={label}>
        {label}
      </label>
      <input
        type="number"
        step="any"
        value={value ?? ''}
        onChange={(e) =>
          onChange(e.target.value === '' ? null : parseFloat(e.target.value))
        }
        placeholder="--"
        className="flex-1 px-2 py-1 bg-gray-800 border border-gray-700 rounded
                   text-xs font-mono text-gray-200 focus:border-blue-500
                   focus:outline-none focus:ring-1 focus:ring-blue-500/50
                   placeholder:text-gray-600
                   [appearance:textfield]
                   [&::-webkit-outer-spin-button]:appearance-none
                   [&::-webkit-inner-spin-button]:appearance-none"
      />
      {unit && <span className="text-xs text-gray-500 w-16 shrink-0">{unit}</span>}
    </div>
  );
}

function ComponentSelector({
  label,
  value,
  options,
  onChange,
  color,
}: {
  label: string;
  value: string | null;
  options: { id: string; name: string; detail: string }[];
  onChange: (id: string | null) => void;
  color: string;
}) {
  const isAssigned = value !== null;
  return (
    <div className="bg-gray-900/50 rounded-lg p-3 border border-gray-800">
      <div className="flex items-center gap-2 mb-2">
        {isAssigned ? (
          <CheckCircle2 size={14} className={`text-${color}-400`} />
        ) : (
          <XCircle size={14} className="text-gray-600" />
        )}
        <span className="text-xs font-semibold text-gray-300">{label}</span>
      </div>
      <select
        value={value ?? ''}
        onChange={(e) => onChange(e.target.value || null)}
        className={clsx(
          'w-full px-2 py-1.5 bg-gray-800 border rounded text-xs text-gray-200',
          'focus:outline-none focus:ring-1',
          isAssigned
            ? `border-${color}-700 focus:border-${color}-500 focus:ring-${color}-500/50`
            : 'border-gray-700 focus:border-gray-500 focus:ring-gray-500/50'
        )}
      >
        <option value="">-- Not Assigned --</option>
        {options.map((opt) => (
          <option key={opt.id} value={opt.id}>
            {opt.name} ({opt.detail})
          </option>
        ))}
      </select>
    </div>
  );
}

function SectionDivider({ title }: { title: string }) {
  return (
    <div className="flex items-center gap-2 mt-5 mb-2">
      <span className="text-xs font-bold text-gray-400 uppercase tracking-wider">{title}</span>
      <div className="flex-1 border-b border-gray-800" />
    </div>
  );
}

// ---------------------------------------------------------------------------
// Turbine Schematic SVG
// ---------------------------------------------------------------------------

function TurbineSchematic({
  model,
  towerHeight,
  bladeLength,
}: {
  model: TurbineModel;
  towerHeight: number | null;
  bladeLength: number | null;
}) {
  const svgWidth = 500;
  const svgHeight = 550;
  const padTop = 40;
  const padBottom = 50;
  const drawH = svgHeight - padTop - padBottom;

  const tH = towerHeight ?? 87.6;
  const bL = bladeLength ?? 61.5;
  const overhang = model.overhang ?? 5;
  const shaftTilt = model.shaft_tilt ?? 5;
  const precone = model.precone ?? 2.5;

  const totalHeight = tH + bL * 1.1;
  const scale = drawH / totalHeight;

  const towerBaseY = svgHeight - padBottom;
  const towerTopY = towerBaseY - tH * scale;
  const towerBaseDia = 6 * scale;
  const towerTopDia = 3.87 * scale;
  const centerX = svgWidth * 0.45;

  // Nacelle
  const nacelleW = 12 * scale;
  const nacelleH = 3.5 * scale;
  const nacelleX = centerX + overhang * scale - nacelleW * 0.3;
  const nacelleY = towerTopY - nacelleH;

  // Hub center
  const hubX = nacelleX + nacelleW;
  const hubY = nacelleY + nacelleH * 0.5;
  const hubR = 2 * scale;

  // Blades (3 blades, 120 degrees apart)
  const bladeLenPx = bL * scale;
  const tiltRad = (shaftTilt * Math.PI) / 180;
  const preconeRad = (precone * Math.PI) / 180;

  const bladeAngles = [90, 210, 330]; // degrees from vertical
  const bladePaths = bladeAngles.map((angle) => {
    const rad = (angle * Math.PI) / 180 - tiltRad;
    const tipX = hubX + Math.sin(rad + preconeRad) * bladeLenPx;
    const tipY = hubY - Math.cos(rad + preconeRad) * bladeLenPx;
    // blade width tapers
    const rootW = 4 * scale;
    const tipW = 1 * scale;
    const perpX = Math.cos(rad);
    const perpY = Math.sin(rad);
    return {
      path: [
        `M${hubX + perpX * rootW},${hubY + perpY * rootW}`,
        `L${tipX + perpX * tipW},${tipY + perpY * tipW}`,
        `L${tipX - perpX * tipW},${tipY - perpY * tipW}`,
        `L${hubX - perpX * rootW},${hubY - perpY * rootW}`,
        'Z',
      ].join(' '),
      tipX,
      tipY,
    };
  });

  // Rotor diameter dimension
  const rotorDiameter = bL * 2;
  const topBladeTipY = Math.min(...bladePaths.map((b) => b.tipY));

  return (
    <svg
      viewBox={`0 0 ${svgWidth} ${svgHeight}`}
      className="w-full h-full"
      style={{ maxHeight: '100%' }}
    >
      <defs>
        <pattern id="turbGrid" width="25" height="25" patternUnits="userSpaceOnUse">
          <path d="M 25 0 L 0 0 0 25" fill="none" stroke="rgba(75,85,99,0.12)" strokeWidth="0.5" />
        </pattern>
        <linearGradient id="towerBodyGrad" x1="0" y1="0" x2="1" y2="0">
          <stop offset="0%" stopColor="rgba(59,130,246,0.15)" />
          <stop offset="50%" stopColor="rgba(59,130,246,0.3)" />
          <stop offset="100%" stopColor="rgba(59,130,246,0.15)" />
        </linearGradient>
        <linearGradient id="nacelleGrad" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0%" stopColor="rgba(156,163,175,0.4)" />
          <stop offset="100%" stopColor="rgba(75,85,99,0.4)" />
        </linearGradient>
        <linearGradient id="bladeGrad" x1="0" y1="0" x2="1" y2="1">
          <stop offset="0%" stopColor="rgba(16,185,129,0.2)" />
          <stop offset="100%" stopColor="rgba(16,185,129,0.35)" />
        </linearGradient>
      </defs>
      <rect width={svgWidth} height={svgHeight} fill="url(#turbGrid)" />

      {/* Ground */}
      <line x1={0} y1={towerBaseY} x2={svgWidth} y2={towerBaseY} stroke="#6b7280" strokeWidth={2} />
      <text x={svgWidth - 8} y={towerBaseY + 14} textAnchor="end" fontSize={9} className="fill-gray-500">
        Ground
      </text>

      {/* Tower */}
      <path
        d={`M${centerX - towerBaseDia / 2},${towerBaseY}
            L${centerX - towerTopDia / 2},${towerTopY}
            L${centerX + towerTopDia / 2},${towerTopY}
            L${centerX + towerBaseDia / 2},${towerBaseY} Z`}
        fill="url(#towerBodyGrad)"
        stroke="#3b82f6"
        strokeWidth={1.5}
      />

      {/* Nacelle */}
      <rect
        x={nacelleX}
        y={nacelleY}
        width={nacelleW}
        height={nacelleH}
        rx={2}
        fill="url(#nacelleGrad)"
        stroke="#9ca3af"
        strokeWidth={1}
      />

      {/* Blades */}
      {bladePaths.map((bp, i) => (
        <path
          key={i}
          d={bp.path}
          fill="url(#bladeGrad)"
          stroke="#10b981"
          strokeWidth={1}
        />
      ))}

      {/* Hub */}
      <circle cx={hubX} cy={hubY} r={hubR} fill="#374151" stroke="#9ca3af" strokeWidth={1.5} />

      {/* ---- Dimension lines ---- */}

      {/* Tower height */}
      <line
        x1={centerX - towerBaseDia / 2 - 20}
        y1={towerBaseY}
        x2={centerX - towerBaseDia / 2 - 20}
        y2={towerTopY}
        stroke="#6b7280"
        strokeWidth={0.5}
      />
      <line
        x1={centerX - towerBaseDia / 2 - 25}
        y1={towerBaseY}
        x2={centerX - towerBaseDia / 2 - 15}
        y2={towerBaseY}
        stroke="#6b7280"
        strokeWidth={0.5}
      />
      <line
        x1={centerX - towerBaseDia / 2 - 25}
        y1={towerTopY}
        x2={centerX - towerBaseDia / 2 - 15}
        y2={towerTopY}
        stroke="#6b7280"
        strokeWidth={0.5}
      />
      <text
        x={centerX - towerBaseDia / 2 - 28}
        y={(towerBaseY + towerTopY) / 2}
        textAnchor="end"
        fontSize={9}
        className="fill-gray-400"
        transform={`rotate(-90, ${centerX - towerBaseDia / 2 - 28}, ${(towerBaseY + towerTopY) / 2})`}
      >
        Hub Height: {tH.toFixed(1)} m
      </text>

      {/* Rotor diameter - horizontal line through hub */}
      {bladePaths.length > 0 && (
        <>
          <line
            x1={hubX}
            y1={topBladeTipY - 5}
            x2={hubX}
            y2={topBladeTipY - 15}
            stroke="#10b981"
            strokeWidth={0.5}
          />
          <text
            x={hubX}
            y={topBladeTipY - 20}
            textAnchor="middle"
            fontSize={9}
            className="fill-green-400"
          >
            Rotor: {rotorDiameter.toFixed(1)} m
          </text>
        </>
      )}

      {/* Overhang dimension */}
      {overhang > 0 && (
        <>
          <line
            x1={centerX}
            y1={towerTopY - 5}
            x2={hubX}
            y2={towerTopY - 5}
            stroke="#f59e0b"
            strokeWidth={0.5}
            strokeDasharray="2,2"
          />
          <text
            x={(centerX + hubX) / 2}
            y={towerTopY - 10}
            textAnchor="middle"
            fontSize={8}
            className="fill-amber-400"
          >
            Overhang: {overhang.toFixed(1)} m
          </text>
        </>
      )}

      {/* Shaft tilt arc indicator */}
      {shaftTilt > 0 && (
        <text
          x={hubX + 10}
          y={hubY - 5}
          fontSize={8}
          className="fill-gray-500"
        >
          Tilt: {shaftTilt.toFixed(1)} deg
        </text>
      )}

      {/* Precone label */}
      {precone > 0 && (
        <text
          x={hubX + bladeLenPx * 0.3}
          y={hubY - bladeLenPx * 0.3 - 5}
          fontSize={8}
          className="fill-gray-500"
        >
          Precone: {precone.toFixed(1)} deg
        </text>
      )}
    </svg>
  );
}

// ---------------------------------------------------------------------------
// Main TurbineAssembly
// ---------------------------------------------------------------------------

export default function TurbineAssembly() {
  const { projectId } = useParams<{ projectId: string }>();

  const [models, setModels] = useState<TurbineModel[]>([]);
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [model, setModel] = useState<TurbineModel | null>(null);
  const [originalModel, setOriginalModel] = useState<TurbineModel | null>(null);

  const [towers, setTowers] = useState<Tower[]>([]);
  const [blades, setBlades] = useState<Blade[]>([]);
  const [controllers, setControllers] = useState<Controller[]>([]);

  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [validationResults, setValidationResults] = useState<
    { field: string; ok: boolean; msg: string }[] | null
  >(null);
  const [dropdownOpen, setDropdownOpen] = useState(false);

  const isDirty = useMemo(
    () => JSON.stringify(model) !== JSON.stringify(originalModel),
    [model, originalModel]
  );

  // Selected component details
  const selectedTower = useMemo(
    () => towers.find((t) => t.id === model?.tower_id) ?? null,
    [towers, model?.tower_id]
  );
  const selectedBlade = useMemo(
    () => blades.find((b) => b.id === model?.blade_id) ?? null,
    [blades, model?.blade_id]
  );

  // ---------------------------------------------------------------------------
  // Fetching
  // ---------------------------------------------------------------------------

  const fetchAll = useCallback(async () => {
    if (!projectId) return;
    try {
      setLoading(true);
      const [modelsResp, towersResp, bladesResp, ctrlResp] = await Promise.all([
        axios.get<TurbineModel[]>(`${API_BASE}/projects/${projectId}/turbine-models`),
        axios.get<Tower[]>(`${API_BASE}/projects/${projectId}/towers`),
        axios.get<Blade[]>(`${API_BASE}/projects/${projectId}/blades`),
        axios.get<Controller[]>(`${API_BASE}/projects/${projectId}/controllers`),
      ]);
      setModels(modelsResp.data);
      setTowers(towersResp.data);
      setBlades(bladesResp.data);
      setControllers(ctrlResp.data);
      if (modelsResp.data.length > 0 && !selectedId) {
        setSelectedId(modelsResp.data[0].id);
      }
    } catch (err: any) {
      setError(err.message || 'Failed to load data');
    } finally {
      setLoading(false);
    }
  }, [projectId, selectedId]);

  const fetchModel = useCallback(
    async (modelId: string) => {
      if (!projectId) return;
      try {
        const resp = await axios.get<TurbineModel>(
          `${API_BASE}/projects/${projectId}/turbine-models/${modelId}`
        );
        setModel(resp.data);
        setOriginalModel(JSON.parse(JSON.stringify(resp.data)));
        setError(null);
      } catch (err: any) {
        setError(err.message || 'Failed to load turbine model');
      }
    },
    [projectId]
  );

  useEffect(() => {
    fetchAll();
  }, [fetchAll]);

  useEffect(() => {
    if (selectedId) fetchModel(selectedId);
  }, [selectedId, fetchModel]);

  // ---------------------------------------------------------------------------
  // Actions
  // ---------------------------------------------------------------------------

  const handleSave = useCallback(async () => {
    if (!projectId || !model) return;
    try {
      setSaving(true);
      if (model.id) {
        await axios.put(
          `${API_BASE}/projects/${projectId}/turbine-models/${model.id}`,
          model
        );
      } else {
        const resp = await axios.post<TurbineModel>(
          `${API_BASE}/projects/${projectId}/turbine-models`,
          model
        );
        setSelectedId(resp.data.id);
      }
      setOriginalModel(JSON.parse(JSON.stringify(model)));
      await fetchAll();
      setError(null);
    } catch (err: any) {
      setError(err.message || 'Failed to save turbine model');
    } finally {
      setSaving(false);
    }
  }, [projectId, model, fetchAll]);

  const handleRevert = useCallback(() => {
    if (originalModel) setModel(JSON.parse(JSON.stringify(originalModel)));
    setValidationResults(null);
  }, [originalModel]);

  const handleNewModel = useCallback(() => {
    const newModel: TurbineModel = {
      id: '',
      project_id: projectId || '',
      name: 'New Turbine Model',
      version: 1,
      tower_id: null,
      blade_id: null,
      controller_id: null,
      gearbox_ratio: 97.0,
      generator_inertia: 534.116,
      drivetrain_stiffness: 867637000,
      drivetrain_damping: 6215000,
      hub_mass: 56780,
      hub_inertia: 115926,
      nacelle_mass: 240000,
      nacelle_inertia: 2607890,
      overhang: -5.0191,
      shaft_tilt: 5.0,
      precone: -2.5,
      rotor_speed_rated: 12.1,
      dof_flags: {
        FlapDOF1: true,
        FlapDOF2: true,
        EdgeDOF1: true,
        TwrFADOF1: true,
        TwrFADOF2: true,
        TwrSSDOF1: true,
        TwrSSDOF2: true,
        DrTrDOF: true,
        GenDOF: true,
        YawDOF: false,
        PtfmSgDOF: false,
        PtfmSwDOF: false,
        PtfmHvDOF: false,
      },
      is_active: true,
      created_at: new Date().toISOString(),
    };
    setModel(newModel);
    setOriginalModel(null);
    setSelectedId(null);
    setValidationResults(null);
  }, [projectId]);

  const updateField = useCallback(
    <K extends keyof TurbineModel>(key: K, value: TurbineModel[K]) => {
      if (!model) return;
      setModel({ ...model, [key]: value });
    },
    [model]
  );

  const updateDofFlag = useCallback(
    (key: string, value: boolean) => {
      if (!model) return;
      setModel({
        ...model,
        dof_flags: { ...(model.dof_flags ?? {}), [key]: value },
      });
    },
    [model]
  );

  const handleValidate = useCallback(() => {
    if (!model) return;
    const results: { field: string; ok: boolean; msg: string }[] = [];

    results.push({
      field: 'Tower',
      ok: model.tower_id !== null,
      msg: model.tower_id ? 'Assigned' : 'No tower assigned',
    });
    results.push({
      field: 'Blade',
      ok: model.blade_id !== null,
      msg: model.blade_id ? 'Assigned' : 'No blade assigned',
    });
    results.push({
      field: 'Controller',
      ok: model.controller_id !== null,
      msg: model.controller_id ? 'Assigned' : 'No controller assigned',
    });
    results.push({
      field: 'Gearbox Ratio',
      ok: model.gearbox_ratio !== null && model.gearbox_ratio > 0,
      msg: model.gearbox_ratio ? `${model.gearbox_ratio}` : 'Not set',
    });
    results.push({
      field: 'Generator Inertia',
      ok: model.generator_inertia !== null && model.generator_inertia > 0,
      msg: model.generator_inertia ? `${model.generator_inertia} kg-m^2` : 'Not set',
    });
    results.push({
      field: 'Hub Mass',
      ok: model.hub_mass !== null && model.hub_mass > 0,
      msg: model.hub_mass ? `${model.hub_mass} kg` : 'Not set',
    });
    results.push({
      field: 'Nacelle Mass',
      ok: model.nacelle_mass !== null && model.nacelle_mass > 0,
      msg: model.nacelle_mass ? `${model.nacelle_mass} kg` : 'Not set',
    });
    results.push({
      field: 'Rotor Speed',
      ok: model.rotor_speed_rated !== null && model.rotor_speed_rated > 0,
      msg: model.rotor_speed_rated ? `${model.rotor_speed_rated} rpm` : 'Not set',
    });

    setValidationResults(results);
  }, [model]);

  // ---------------------------------------------------------------------------
  // Render
  // ---------------------------------------------------------------------------

  if (loading && !model) {
    return (
      <div className="flex items-center justify-center h-full bg-gray-950 text-gray-400">
        <Loader2 className="animate-spin mr-2" size={20} />
        Loading turbine model data...
      </div>
    );
  }

  return (
    <div className="flex flex-col h-full bg-gray-950 text-gray-200">
      {/* Top bar */}
      <div className="flex items-center gap-3 px-4 py-2 bg-gray-900 border-b border-gray-800 shrink-0">
        <Zap size={16} className="text-amber-400" />
        <span className="text-sm font-semibold text-gray-200">Turbine Assembly</span>

        {/* Model selector */}
        <div className="relative ml-4">
          <button
            onClick={() => setDropdownOpen(!dropdownOpen)}
            className="flex items-center gap-2 px-3 py-1.5 bg-gray-800 border border-gray-700
                       rounded text-xs hover:bg-gray-750 transition-colors min-w-[200px]"
          >
            <span className="truncate">{model?.name || 'Select Model'}</span>
            <ChevronDown size={14} className="ml-auto text-gray-500" />
          </button>
          {dropdownOpen && (
            <>
              <div className="fixed inset-0 z-10" onClick={() => setDropdownOpen(false)} />
              <div className="absolute top-full left-0 mt-1 w-64 bg-gray-800 border border-gray-700
                              rounded-lg shadow-xl z-20 max-h-60 overflow-auto">
                {models.map((m) => (
                  <button
                    key={m.id}
                    onClick={() => {
                      setSelectedId(m.id);
                      setDropdownOpen(false);
                    }}
                    className={clsx(
                      'w-full text-left px-3 py-2 text-xs hover:bg-gray-700 transition-colors',
                      m.id === selectedId && 'bg-amber-900/40 text-amber-300'
                    )}
                  >
                    <div className="font-medium">{m.name}</div>
                    <div className="text-gray-500">v{m.version}</div>
                  </button>
                ))}
                {models.length === 0 && (
                  <div className="px-3 py-4 text-xs text-gray-500 text-center">
                    No turbine models yet
                  </div>
                )}
              </div>
            </>
          )}
        </div>

        <button
          onClick={handleNewModel}
          className="flex items-center gap-1 px-3 py-1.5 text-xs bg-gray-800
                     border border-gray-700 rounded hover:bg-gray-750 transition-colors"
        >
          <Plus size={14} />
          New Model
        </button>

        {model && (
          <span className="ml-2 px-2 py-0.5 bg-amber-900/40 text-amber-300 text-xs rounded-full font-medium">
            v{model.version}
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
            disabled={saving || !model}
            className="flex items-center gap-1 px-4 py-1.5 text-xs bg-amber-600
                       text-white rounded hover:bg-amber-700
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
      {model ? (
        <div className="flex-1 flex overflow-hidden">
          {/* LEFT PANEL: Component selection (~40%) */}
          <div className="w-[40%] flex flex-col border-r border-gray-800 overflow-hidden">
            <div className="p-3 space-y-3 overflow-auto flex-1">
              {/* Model name */}
              <div className="flex items-center gap-2">
                <label className="text-xs text-gray-400 w-20 shrink-0">Name</label>
                <input
                  type="text"
                  value={model.name}
                  onChange={(e) => updateField('name', e.target.value)}
                  className="flex-1 px-2 py-1.5 bg-gray-800 border border-gray-700 rounded
                             text-sm text-gray-200 focus:border-amber-500
                             focus:outline-none focus:ring-1 focus:ring-amber-500/50"
                />
              </div>

              {/* Component selectors */}
              <SectionDivider title="Components" />
              <div className="space-y-2">
                <ComponentSelector
                  label="Tower"
                  value={model.tower_id}
                  options={towers.map((t) => ({
                    id: t.id,
                    name: t.name,
                    detail: `v${t.version}, ${t.tower_height}m`,
                  }))}
                  onChange={(id) => updateField('tower_id', id)}
                  color="blue"
                />
                <ComponentSelector
                  label="Blade"
                  value={model.blade_id}
                  options={blades.map((b) => ({
                    id: b.id,
                    name: b.name,
                    detail: `v${b.version}, ${b.blade_length}m`,
                  }))}
                  onChange={(id) => updateField('blade_id', id)}
                  color="green"
                />
                <ComponentSelector
                  label="Controller"
                  value={model.controller_id}
                  options={controllers.map((c) => ({
                    id: c.id,
                    name: c.name,
                    detail: `v${c.version}, ${c.controller_type}`,
                  }))}
                  onChange={(id) => updateField('controller_id', id)}
                  color="purple"
                />
              </div>

              {/* Drivetrain parameters */}
              <SectionDivider title="Drivetrain" />
              <div className="bg-gray-900/50 rounded-lg p-3 space-y-2 border border-gray-800">
                <ParamInput
                  label="Gearbox Ratio"
                  value={model.gearbox_ratio}
                  unit="-"
                  onChange={(v) => updateField('gearbox_ratio', v)}
                />
                <ParamInput
                  label="Generator Inertia"
                  value={model.generator_inertia}
                  unit="kg-m^2"
                  onChange={(v) => updateField('generator_inertia', v)}
                />
                <ParamInput
                  label="Drivetrain Stiffness"
                  value={model.drivetrain_stiffness}
                  unit="N-m/rad"
                  onChange={(v) => updateField('drivetrain_stiffness', v)}
                />
                <ParamInput
                  label="Drivetrain Damping"
                  value={model.drivetrain_damping}
                  unit="N-m/(rad/s)"
                  onChange={(v) => updateField('drivetrain_damping', v)}
                />
              </div>

              {/* Hub / Nacelle */}
              <SectionDivider title="Hub / Nacelle" />
              <div className="bg-gray-900/50 rounded-lg p-3 space-y-2 border border-gray-800">
                <ParamInput
                  label="Hub Mass"
                  value={model.hub_mass}
                  unit="kg"
                  onChange={(v) => updateField('hub_mass', v)}
                />
                <ParamInput
                  label="Hub Inertia"
                  value={model.hub_inertia}
                  unit="kg-m^2"
                  onChange={(v) => updateField('hub_inertia', v)}
                />
                <ParamInput
                  label="Nacelle Mass"
                  value={model.nacelle_mass}
                  unit="kg"
                  onChange={(v) => updateField('nacelle_mass', v)}
                />
                <ParamInput
                  label="Nacelle Inertia"
                  value={model.nacelle_inertia}
                  unit="kg-m^2"
                  onChange={(v) => updateField('nacelle_inertia', v)}
                />
                <ParamInput
                  label="Overhang"
                  value={model.overhang}
                  unit="m"
                  onChange={(v) => updateField('overhang', v)}
                />
                <ParamInput
                  label="Shaft Tilt"
                  value={model.shaft_tilt}
                  unit="deg"
                  onChange={(v) => updateField('shaft_tilt', v)}
                />
                <ParamInput
                  label="Precone"
                  value={model.precone}
                  unit="deg"
                  onChange={(v) => updateField('precone', v)}
                />
                <ParamInput
                  label="Rated Rotor Speed"
                  value={model.rotor_speed_rated}
                  unit="rpm"
                  onChange={(v) => updateField('rotor_speed_rated', v)}
                />
              </div>

              {/* DOF Flags */}
              <SectionDivider title="Degrees of Freedom" />
              <div className="bg-gray-900/50 rounded-lg p-3 border border-gray-800">
                <div className="grid grid-cols-2 gap-x-4 gap-y-1.5">
                  {DOF_FLAGS.map((dof) => (
                    <label
                      key={dof.key}
                      className="flex items-center gap-2 text-xs text-gray-300 cursor-pointer
                                 hover:text-gray-100 transition-colors"
                    >
                      <input
                        type="checkbox"
                        checked={model.dof_flags?.[dof.key] ?? false}
                        onChange={(e) => updateDofFlag(dof.key, e.target.checked)}
                        className="accent-amber-500 w-3.5 h-3.5"
                      />
                      <span className="font-mono text-[11px]">{dof.label}</span>
                    </label>
                  ))}
                </div>
              </div>
            </div>
          </div>

          {/* RIGHT PANEL: Visualization (~60%) */}
          <div className="w-[60%] flex flex-col overflow-hidden">
            <div className="px-3 py-2 bg-gray-900/50 border-b border-gray-800 shrink-0 flex items-center">
              <span className="text-xs font-semibold text-gray-300">Turbine Schematic</span>
              <button
                onClick={handleValidate}
                className="ml-auto flex items-center gap-1 px-3 py-1.5 text-xs
                           bg-gray-800 border border-gray-700 rounded
                           hover:bg-gray-750 transition-colors"
              >
                <ShieldCheck size={14} />
                Validate All
              </button>
            </div>

            <div className="flex-1 overflow-auto">
              {/* Schematic */}
              <div className="p-4 flex items-center justify-center" style={{ minHeight: '450px' }}>
                <TurbineSchematic
                  model={model}
                  towerHeight={selectedTower?.tower_height ?? null}
                  bladeLength={selectedBlade?.blade_length ?? null}
                />
              </div>

              {/* Component status badges */}
              <div className="px-4 pb-4">
                <div className="grid grid-cols-3 gap-3">
                  {/* Tower status */}
                  <div
                    className={clsx(
                      'rounded-lg p-3 border',
                      model.tower_id
                        ? 'border-blue-800 bg-blue-900/20'
                        : 'border-gray-800 bg-gray-900/50'
                    )}
                  >
                    <div className="flex items-center gap-2 mb-1">
                      {model.tower_id ? (
                        <CheckCircle2 size={14} className="text-blue-400" />
                      ) : (
                        <XCircle size={14} className="text-gray-600" />
                      )}
                      <span className="text-xs font-semibold text-gray-300">Tower</span>
                    </div>
                    <div className="text-xs text-gray-500">
                      {selectedTower
                        ? `${selectedTower.name} (${selectedTower.tower_height}m)`
                        : 'Not assigned'}
                    </div>
                  </div>

                  {/* Blade status */}
                  <div
                    className={clsx(
                      'rounded-lg p-3 border',
                      model.blade_id
                        ? 'border-green-800 bg-green-900/20'
                        : 'border-gray-800 bg-gray-900/50'
                    )}
                  >
                    <div className="flex items-center gap-2 mb-1">
                      {model.blade_id ? (
                        <CheckCircle2 size={14} className="text-green-400" />
                      ) : (
                        <XCircle size={14} className="text-gray-600" />
                      )}
                      <span className="text-xs font-semibold text-gray-300">Blade</span>
                    </div>
                    <div className="text-xs text-gray-500">
                      {selectedBlade
                        ? `${selectedBlade.name} (${selectedBlade.blade_length}m)`
                        : 'Not assigned'}
                    </div>
                  </div>

                  {/* Controller status */}
                  <div
                    className={clsx(
                      'rounded-lg p-3 border',
                      model.controller_id
                        ? 'border-purple-800 bg-purple-900/20'
                        : 'border-gray-800 bg-gray-900/50'
                    )}
                  >
                    <div className="flex items-center gap-2 mb-1">
                      {model.controller_id ? (
                        <CheckCircle2 size={14} className="text-purple-400" />
                      ) : (
                        <XCircle size={14} className="text-gray-600" />
                      )}
                      <span className="text-xs font-semibold text-gray-300">Controller</span>
                    </div>
                    <div className="text-xs text-gray-500">
                      {controllers.find((c) => c.id === model.controller_id)?.name || 'Not assigned'}
                    </div>
                  </div>
                </div>
              </div>

              {/* Validation results */}
              {validationResults && (
                <div className="px-4 pb-4">
                  <div className="bg-gray-900/50 rounded-lg border border-gray-800 overflow-hidden">
                    <div className="px-3 py-2 bg-gray-800 border-b border-gray-700 flex items-center gap-2">
                      <ShieldCheck size={14} className="text-gray-400" />
                      <span className="text-xs font-semibold text-gray-300">Validation Results</span>
                      <span
                        className={clsx(
                          'ml-auto text-xs font-medium',
                          validationResults.every((r) => r.ok)
                            ? 'text-green-400'
                            : 'text-amber-400'
                        )}
                      >
                        {validationResults.filter((r) => r.ok).length}/{validationResults.length} passed
                      </span>
                    </div>
                    <div className="divide-y divide-gray-800">
                      {validationResults.map((result, i) => (
                        <div
                          key={i}
                          className="flex items-center gap-2 px-3 py-2 text-xs"
                        >
                          {result.ok ? (
                            <CheckCircle2 size={13} className="text-green-400 shrink-0" />
                          ) : (
                            <XCircle size={13} className="text-red-400 shrink-0" />
                          )}
                          <span className="text-gray-300 w-36 shrink-0">{result.field}</span>
                          <span
                            className={clsx(
                              'truncate',
                              result.ok ? 'text-gray-500' : 'text-red-400'
                            )}
                          >
                            {result.msg}
                          </span>
                        </div>
                      ))}
                    </div>
                  </div>
                </div>
              )}
            </div>
          </div>
        </div>
      ) : (
        <div className="flex-1 flex items-center justify-center">
          <div className="text-center text-gray-500">
            <Zap size={64} className="mx-auto mb-4 opacity-20" />
            <p className="text-lg font-medium">No turbine model selected</p>
            <p className="text-sm mt-2">Select a model or create a new one</p>
            <button
              onClick={handleNewModel}
              className="mt-4 px-4 py-2 bg-amber-600 text-white rounded
                         hover:bg-amber-700 transition-colors text-sm"
            >
              Create New Model
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
