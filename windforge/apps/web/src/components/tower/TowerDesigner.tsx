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
  Box,
} from 'lucide-react';
import clsx from 'clsx';
import DataGrid, { DataGridColumn } from '@/components/common/DataGrid';
import FilePreview from '@/components/common/FilePreview';
import PlotPanel from '@/components/common/PlotPanel';

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface TowerStation {
  frac: number;
  mass_den: number;
  fa_stiff: number;
  ss_stiff: number;
  outer_diameter: number;
  wall_thickness: number;
}

interface Tower {
  id: string;
  project_id: string;
  name: string;
  version: number;
  tower_height: number;
  tower_base_height: number;
  tower_fa_damping_1: number;
  tower_fa_damping_2: number;
  tower_ss_damping_1: number;
  tower_ss_damping_2: number;
  stations: TowerStation[] | null;
  fa_mode_1_coeffs: number[] | null;
  fa_mode_2_coeffs: number[] | null;
  ss_mode_1_coeffs: number[] | null;
  ss_mode_2_coeffs: number[] | null;
  is_active: boolean;
  created_at: string;
}

const API_BASE = '/api/v1';

const STATION_COLUMNS: DataGridColumn[] = [
  { name: 'Frac', key: 'frac', type: 'number', min: 0, max: 1, step: 0.01, unit: '-', precision: 4 },
  { name: 'Mass Den', key: 'mass_den', type: 'number', min: 0, step: 1, unit: 'kg/m', precision: 2 },
  { name: 'FA Stiff', key: 'fa_stiff', type: 'number', min: 0, step: 1e6, unit: 'N-m^2', precision: 2 },
  { name: 'SS Stiff', key: 'ss_stiff', type: 'number', min: 0, step: 1e6, unit: 'N-m^2', precision: 2 },
  { name: 'OD', key: 'outer_diameter', type: 'number', min: 0, step: 0.1, unit: 'm', precision: 4 },
  { name: 'Wall Thick', key: 'wall_thickness', type: 'number', min: 0, step: 0.001, unit: 'm', precision: 4 },
];

const DEFAULT_STATION: TowerStation = {
  frac: 0,
  mass_den: 0,
  fa_stiff: 0,
  ss_stiff: 0,
  outer_diameter: 0,
  wall_thickness: 0,
};

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
  value: number;
  unit?: string;
  onChange: (v: number) => void;
}) {
  return (
    <div className="flex items-center gap-2">
      <label className="text-xs text-gray-400 w-32 shrink-0">{label}</label>
      <input
        type="number"
        step="any"
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value) || 0)}
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

function ModeShapeInputs({
  label,
  coeffs,
  onChange,
}: {
  label: string;
  coeffs: number[];
  onChange: (c: number[]) => void;
}) {
  return (
    <div className="space-y-1">
      <div className="text-xs font-semibold text-gray-400">{label}</div>
      <div className="grid grid-cols-5 gap-1">
        {coeffs.map((c, i) => (
          <input
            key={i}
            type="number"
            step="any"
            value={c}
            onChange={(e) => {
              const next = [...coeffs];
              next[i] = parseFloat(e.target.value) || 0;
              onChange(next);
            }}
            title={`x^${i + 2}`}
            className="px-1 py-1 bg-gray-800 border border-gray-700 rounded
                       text-xs font-mono text-gray-200 text-center
                       focus:border-blue-500 focus:outline-none
                       [appearance:textfield]
                       [&::-webkit-outer-spin-button]:appearance-none
                       [&::-webkit-inner-spin-button]:appearance-none"
          />
        ))}
      </div>
      <div className="grid grid-cols-5 gap-1">
        {coeffs.map((_, i) => (
          <span key={i} className="text-[9px] text-gray-600 text-center">
            x^{i + 2}
          </span>
        ))}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Tower SVG Visualization
// ---------------------------------------------------------------------------

function TowerVisualization({
  stations,
  towerHeight,
  towerBaseHeight,
}: {
  stations: TowerStation[];
  towerHeight: number;
  towerBaseHeight: number;
}) {
  const svgWidth = 280;
  const svgHeight = 500;
  const padTop = 30;
  const padBottom = 50;
  const padX = 60;
  const drawHeight = svgHeight - padTop - padBottom;

  const maxOD = useMemo(
    () => Math.max(...stations.map((s) => s.outer_diameter), 1),
    [stations]
  );
  const maxDrawWidth = svgWidth - padX * 2;
  const scale = maxDrawWidth / maxOD;
  const centerX = svgWidth / 2;

  const sortedStations = useMemo(
    () => [...stations].sort((a, b) => a.frac - b.frac),
    [stations]
  );

  // Build tower outline points
  const leftPoints: [number, number][] = [];
  const rightPoints: [number, number][] = [];

  sortedStations.forEach((st) => {
    const y = padTop + drawHeight * (1 - st.frac);
    const halfW = (st.outer_diameter * scale) / 2;
    leftPoints.push([centerX - halfW, y]);
    rightPoints.push([centerX + halfW, y]);
  });

  // Polygon outline
  const outlinePoints = [
    ...leftPoints,
    ...rightPoints.reverse(),
  ];

  const outlinePath = outlinePoints.length > 0
    ? outlinePoints.map((p, i) => `${i === 0 ? 'M' : 'L'}${p[0]},${p[1]}`).join(' ') + ' Z'
    : '';

  // Height labels
  const heightLabels = useMemo(() => {
    if (sortedStations.length === 0) return [];
    const labels: { frac: number; height: number }[] = [];
    const numLabels = Math.min(sortedStations.length, 6);
    const step = Math.max(1, Math.floor((sortedStations.length - 1) / (numLabels - 1)));
    for (let i = 0; i < sortedStations.length; i += step) {
      labels.push({
        frac: sortedStations[i].frac,
        height: towerBaseHeight + sortedStations[i].frac * towerHeight,
      });
    }
    // Always include last
    const last = sortedStations[sortedStations.length - 1];
    if (labels[labels.length - 1]?.frac !== last.frac) {
      labels.push({
        frac: last.frac,
        height: towerBaseHeight + last.frac * towerHeight,
      });
    }
    return labels;
  }, [sortedStations, towerHeight, towerBaseHeight]);

  return (
    <svg
      viewBox={`0 0 ${svgWidth} ${svgHeight}`}
      className="w-full h-full"
      style={{ maxHeight: '100%' }}
    >
      {/* Background grid */}
      <defs>
        <pattern id="towerGrid" width="20" height="20" patternUnits="userSpaceOnUse">
          <path d="M 20 0 L 0 0 0 20" fill="none" stroke="rgba(75,85,99,0.15)" strokeWidth="0.5" />
        </pattern>
      </defs>
      <rect width={svgWidth} height={svgHeight} fill="url(#towerGrid)" />

      {/* Ground line */}
      <line
        x1={0}
        y1={padTop + drawHeight}
        x2={svgWidth}
        y2={padTop + drawHeight}
        stroke="#6b7280"
        strokeWidth={2}
      />
      <text
        x={svgWidth - 5}
        y={padTop + drawHeight + 14}
        textAnchor="end"
        className="fill-gray-500"
        fontSize={9}
      >
        Ground ({towerBaseHeight.toFixed(1)} m)
      </text>

      {/* Tower body */}
      {outlinePath && (
        <>
          <path
            d={outlinePath}
            fill="url(#towerGradient)"
            stroke="#3b82f6"
            strokeWidth={1.5}
          />
          <defs>
            <linearGradient id="towerGradient" x1="0" y1="0" x2="1" y2="0">
              <stop offset="0%" stopColor="rgba(59,130,246,0.15)" />
              <stop offset="50%" stopColor="rgba(59,130,246,0.25)" />
              <stop offset="100%" stopColor="rgba(59,130,246,0.15)" />
            </linearGradient>
          </defs>
        </>
      )}

      {/* Station markers */}
      {sortedStations.map((st, i) => {
        const y = padTop + drawHeight * (1 - st.frac);
        const halfW = (st.outer_diameter * scale) / 2;
        return (
          <g key={i}>
            {/* Station line */}
            <line
              x1={centerX - halfW - 3}
              y1={y}
              x2={centerX + halfW + 3}
              y2={y}
              stroke="rgba(59,130,246,0.4)"
              strokeWidth={0.5}
              strokeDasharray="2,2"
            />
            {/* Left dot */}
            <circle
              cx={centerX - halfW}
              cy={y}
              r={3}
              fill="#3b82f6"
              stroke="#1e3a5f"
              strokeWidth={1}
              className="cursor-pointer hover:r-4"
            />
            {/* Right dot */}
            <circle
              cx={centerX + halfW}
              cy={y}
              r={3}
              fill="#3b82f6"
              stroke="#1e3a5f"
              strokeWidth={1}
              className="cursor-pointer"
            />
            {/* OD dimension for first and last */}
            {(i === 0 || i === sortedStations.length - 1) && (
              <>
                <line
                  x1={centerX - halfW}
                  y1={y + (i === 0 ? 8 : -8)}
                  x2={centerX + halfW}
                  y2={y + (i === 0 ? 8 : -8)}
                  stroke="#9ca3af"
                  strokeWidth={0.5}
                  markerStart="url(#arrowLeft)"
                  markerEnd="url(#arrowRight)"
                />
                <text
                  x={centerX}
                  y={y + (i === 0 ? 18 : -12)}
                  textAnchor="middle"
                  fontSize={8}
                  className="fill-gray-400"
                >
                  {st.outer_diameter.toFixed(2)} m
                </text>
              </>
            )}
          </g>
        );
      })}

      {/* Arrow markers */}
      <defs>
        <marker id="arrowLeft" markerWidth="6" markerHeight="6" refX="6" refY="3" orient="auto">
          <path d="M6,0 L0,3 L6,6" fill="none" stroke="#9ca3af" strokeWidth="1" />
        </marker>
        <marker id="arrowRight" markerWidth="6" markerHeight="6" refX="0" refY="3" orient="auto">
          <path d="M0,0 L6,3 L0,6" fill="none" stroke="#9ca3af" strokeWidth="1" />
        </marker>
      </defs>

      {/* Height labels */}
      {heightLabels.map((hl, i) => {
        const y = padTop + drawHeight * (1 - hl.frac);
        return (
          <g key={i}>
            <line
              x1={8}
              y1={y}
              x2={padX - 15}
              y2={y}
              stroke="rgba(156,163,175,0.3)"
              strokeWidth={0.5}
            />
            <text
              x={5}
              y={y + 3}
              fontSize={8}
              className="fill-gray-500"
              textAnchor="start"
            >
              {hl.height.toFixed(1)}
            </text>
          </g>
        );
      })}

      {/* Height axis label */}
      <text
        x={10}
        y={padTop - 10}
        fontSize={9}
        className="fill-gray-400 font-semibold"
      >
        Height (m)
      </text>
    </svg>
  );
}

// ---------------------------------------------------------------------------
// Main TowerDesigner component
// ---------------------------------------------------------------------------

export default function TowerDesigner() {
  const { projectId } = useParams<{ projectId: string }>();

  // State
  const [towers, setTowers] = useState<Tower[]>([]);
  const [selectedTowerId, setSelectedTowerId] = useState<string | null>(null);
  const [tower, setTower] = useState<Tower | null>(null);
  const [originalTower, setOriginalTower] = useState<Tower | null>(null);
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [filePreview, setFilePreview] = useState<string | null>(null);
  const [filePreviewLoading, setFilePreviewLoading] = useState(false);
  const [rightTab, setRightTab] = useState<'plots' | 'preview'>('plots');
  const [dropdownOpen, setDropdownOpen] = useState(false);

  // Derived
  const stations = tower?.stations ?? [];
  const isDirty = useMemo(
    () => JSON.stringify(tower) !== JSON.stringify(originalTower),
    [tower, originalTower]
  );

  // ---------------------------------------------------------------------------
  // Data fetching
  // ---------------------------------------------------------------------------

  const fetchTowers = useCallback(async () => {
    if (!projectId) return;
    try {
      setLoading(true);
      const resp = await axios.get<Tower[]>(`${API_BASE}/projects/${projectId}/towers`);
      setTowers(resp.data);
      if (resp.data.length > 0 && !selectedTowerId) {
        setSelectedTowerId(resp.data[0].id);
      }
    } catch (err: any) {
      setError(err.message || 'Failed to load towers');
    } finally {
      setLoading(false);
    }
  }, [projectId, selectedTowerId]);

  const fetchTower = useCallback(async (towerId: string) => {
    if (!projectId) return;
    try {
      setLoading(true);
      const resp = await axios.get<Tower>(
        `${API_BASE}/projects/${projectId}/towers/${towerId}`
      );
      setTower(resp.data);
      setOriginalTower(JSON.parse(JSON.stringify(resp.data)));
      setError(null);
    } catch (err: any) {
      setError(err.message || 'Failed to load tower');
    } finally {
      setLoading(false);
    }
  }, [projectId]);

  const fetchPreview = useCallback(async () => {
    if (!projectId || !selectedTowerId) return;
    try {
      setFilePreviewLoading(true);
      const resp = await axios.get<string>(
        `${API_BASE}/projects/${projectId}/towers/${selectedTowerId}/preview`
      );
      setFilePreview(resp.data);
    } catch {
      setFilePreview(null);
    } finally {
      setFilePreviewLoading(false);
    }
  }, [projectId, selectedTowerId]);

  useEffect(() => {
    fetchTowers();
  }, [fetchTowers]);

  useEffect(() => {
    if (selectedTowerId) {
      fetchTower(selectedTowerId);
    }
  }, [selectedTowerId, fetchTower]);

  useEffect(() => {
    if (rightTab === 'preview' && selectedTowerId) {
      fetchPreview();
    }
  }, [rightTab, selectedTowerId, fetchPreview]);

  // ---------------------------------------------------------------------------
  // Actions
  // ---------------------------------------------------------------------------

  const handleSave = useCallback(async () => {
    if (!projectId || !tower) return;
    try {
      setSaving(true);
      if (tower.id) {
        await axios.put(`${API_BASE}/projects/${projectId}/towers/${tower.id}`, tower);
      } else {
        const resp = await axios.post<Tower>(`${API_BASE}/projects/${projectId}/towers`, tower);
        setSelectedTowerId(resp.data.id);
      }
      setOriginalTower(JSON.parse(JSON.stringify(tower)));
      await fetchTowers();
      setError(null);
    } catch (err: any) {
      setError(err.message || 'Failed to save tower');
    } finally {
      setSaving(false);
    }
  }, [projectId, tower, fetchTowers]);

  const handleRevert = useCallback(() => {
    if (originalTower) {
      setTower(JSON.parse(JSON.stringify(originalTower)));
    }
  }, [originalTower]);

  const handleNewTower = useCallback(() => {
    const newTower: Tower = {
      id: '',
      project_id: projectId || '',
      name: 'New Tower',
      version: 1,
      tower_height: 87.6,
      tower_base_height: 10.0,
      tower_fa_damping_1: 1.0,
      tower_fa_damping_2: 1.0,
      tower_ss_damping_1: 1.0,
      tower_ss_damping_2: 1.0,
      stations: [
        { frac: 0.0, mass_den: 5000, fa_stiff: 6.14e11, ss_stiff: 6.14e11, outer_diameter: 6.0, wall_thickness: 0.035 },
        { frac: 0.5, mass_den: 3500, fa_stiff: 2.69e11, ss_stiff: 2.69e11, outer_diameter: 4.935, wall_thickness: 0.025 },
        { frac: 1.0, mass_den: 1800, fa_stiff: 5.49e10, ss_stiff: 5.49e10, outer_diameter: 3.87, wall_thickness: 0.021 },
      ],
      fa_mode_1_coeffs: [0.7004, 2.1963, -5.6202, 6.2275, -2.504],
      fa_mode_2_coeffs: [-17.297, 67.563, -97.851, 69.035, -20.451],
      ss_mode_1_coeffs: [0.7004, 2.1963, -5.6202, 6.2275, -2.504],
      ss_mode_2_coeffs: [-17.297, 67.563, -97.851, 69.035, -20.451],
      is_active: true,
      created_at: new Date().toISOString(),
    };
    setTower(newTower);
    setOriginalTower(null);
    setSelectedTowerId(null);
  }, [projectId]);

  const handleAddStation = useCallback(() => {
    if (!tower) return;
    const current = tower.stations ?? [];
    const lastFrac = current.length > 0 ? current[current.length - 1].frac : 0;
    const newStation: TowerStation = {
      ...DEFAULT_STATION,
      frac: Math.min(lastFrac + 0.1, 1.0),
    };
    setTower({ ...tower, stations: [...current, newStation] });
  }, [tower]);

  const handleDeleteStation = useCallback(
    (index: number) => {
      if (!tower || !tower.stations) return;
      const updated = tower.stations.filter((_, i) => i !== index);
      setTower({ ...tower, stations: updated });
    },
    [tower]
  );

  const handleStationsChange = useCallback(
    (newStations: TowerStation[]) => {
      if (!tower) return;
      setTower({ ...tower, stations: newStations });
    },
    [tower]
  );

  const updateField = useCallback(
    <K extends keyof Tower>(key: K, value: Tower[K]) => {
      if (!tower) return;
      setTower({ ...tower, [key]: value });
    },
    [tower]
  );

  // ---------------------------------------------------------------------------
  // Plot data
  // ---------------------------------------------------------------------------

  const plotData = useMemo(() => {
    const sorted = [...stations].sort((a, b) => a.frac - b.frac);
    const heights = sorted.map(
      (s) => (tower?.tower_base_height ?? 0) + s.frac * (tower?.tower_height ?? 0)
    );

    const massPlot = {
      data: [
        {
          x: heights,
          y: sorted.map((s) => s.mass_den),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'Mass Density',
          line: { color: '#3b82f6', width: 2 },
          marker: { size: 4 },
        },
      ],
      layout: {
        xaxis: { title: { text: 'Height (m)' } },
        yaxis: { title: { text: 'kg/m' } },
      },
      title: 'Mass Density',
    };

    const stiffPlot = {
      data: [
        {
          x: heights,
          y: sorted.map((s) => s.fa_stiff),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'FA Stiffness',
          line: { color: '#10b981', width: 2 },
          marker: { size: 4 },
        },
        {
          x: heights,
          y: sorted.map((s) => s.ss_stiff),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'SS Stiffness',
          line: { color: '#f59e0b', width: 2 },
          marker: { size: 4 },
        },
      ],
      layout: {
        xaxis: { title: { text: 'Height (m)' } },
        yaxis: { title: { text: 'N-m^2' } },
      },
      title: 'Stiffness Distribution',
    };

    const odPlot = {
      data: [
        {
          x: heights,
          y: sorted.map((s) => s.outer_diameter),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'Outer Diameter',
          line: { color: '#8b5cf6', width: 2 },
          marker: { size: 4 },
        },
      ],
      layout: {
        xaxis: { title: { text: 'Height (m)' } },
        yaxis: { title: { text: 'm' } },
      },
      title: 'Outer Diameter',
    };

    return { massPlot, stiffPlot, odPlot };
  }, [stations, tower?.tower_base_height, tower?.tower_height]);

  // ---------------------------------------------------------------------------
  // Render
  // ---------------------------------------------------------------------------

  if (loading && !tower) {
    return (
      <div className="flex items-center justify-center h-full bg-gray-950 text-gray-400">
        <Loader2 className="animate-spin mr-2" size={20} />
        Loading tower data...
      </div>
    );
  }

  return (
    <div className="flex flex-col h-full bg-gray-950 text-gray-200">
      {/* Top bar */}
      <div className="flex items-center gap-3 px-4 py-2 bg-gray-900 border-b border-gray-800 shrink-0">
        <Box size={16} className="text-blue-400" />
        <span className="text-sm font-semibold text-gray-200">Tower Designer</span>

        {/* Tower selector dropdown */}
        <div className="relative ml-4">
          <button
            onClick={() => setDropdownOpen(!dropdownOpen)}
            className="flex items-center gap-2 px-3 py-1.5 bg-gray-800 border border-gray-700
                       rounded text-xs hover:bg-gray-750 transition-colors min-w-[200px]"
          >
            <span className="truncate">
              {tower?.name || 'Select Tower'}
            </span>
            <ChevronDown size={14} className="ml-auto text-gray-500" />
          </button>
          {dropdownOpen && (
            <>
              <div
                className="fixed inset-0 z-10"
                onClick={() => setDropdownOpen(false)}
              />
              <div className="absolute top-full left-0 mt-1 w-64 bg-gray-800 border border-gray-700
                              rounded-lg shadow-xl z-20 max-h-60 overflow-auto">
                {towers.map((t) => (
                  <button
                    key={t.id}
                    onClick={() => {
                      setSelectedTowerId(t.id);
                      setDropdownOpen(false);
                    }}
                    className={clsx(
                      'w-full text-left px-3 py-2 text-xs hover:bg-gray-700 transition-colors',
                      t.id === selectedTowerId && 'bg-blue-900/40 text-blue-300'
                    )}
                  >
                    <div className="font-medium">{t.name}</div>
                    <div className="text-gray-500">v{t.version} &middot; {t.tower_height}m</div>
                  </button>
                ))}
                {towers.length === 0 && (
                  <div className="px-3 py-4 text-xs text-gray-500 text-center">
                    No towers yet
                  </div>
                )}
              </div>
            </>
          )}
        </div>

        <button
          onClick={handleNewTower}
          className="flex items-center gap-1 px-3 py-1.5 text-xs bg-gray-800
                     border border-gray-700 rounded hover:bg-gray-750 transition-colors"
        >
          <Plus size={14} />
          New Tower
        </button>

        {tower && (
          <span className="ml-2 px-2 py-0.5 bg-blue-900/40 text-blue-300 text-xs rounded-full font-medium">
            v{tower.version}
          </span>
        )}

        <div className="ml-auto flex items-center gap-2">
          {isDirty && (
            <span className="text-xs text-amber-400 mr-2">Unsaved changes</span>
          )}
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
            disabled={saving || !tower}
            className="flex items-center gap-1 px-4 py-1.5 text-xs bg-blue-600
                       text-white rounded hover:bg-blue-700
                       disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
          >
            {saving ? (
              <Loader2 size={14} className="animate-spin" />
            ) : (
              <Save size={14} />
            )}
            Save
          </button>
        </div>
      </div>

      {/* Error banner */}
      {error && (
        <div className="flex items-center gap-2 px-4 py-2 bg-red-900/30 border-b border-red-800 text-xs text-red-300">
          <AlertCircle size={14} />
          {error}
          <button
            onClick={() => setError(null)}
            className="ml-auto text-red-400 hover:text-red-200"
          >
            Dismiss
          </button>
        </div>
      )}

      {/* Main content: 3 panels */}
      {tower ? (
        <div className="flex-1 flex overflow-hidden">
          {/* LEFT PANEL: Data table (~40%) */}
          <div className="w-[40%] flex flex-col border-r border-gray-800 overflow-hidden">
            <div className="p-3 space-y-3 overflow-auto flex-1">
              {/* Tower name */}
              <div className="flex items-center gap-2">
                <label className="text-xs text-gray-400 w-20 shrink-0">Name</label>
                <input
                  type="text"
                  value={tower.name}
                  onChange={(e) => updateField('name', e.target.value)}
                  className="flex-1 px-2 py-1.5 bg-gray-800 border border-gray-700 rounded
                             text-sm text-gray-200 focus:border-blue-500
                             focus:outline-none focus:ring-1 focus:ring-blue-500/50"
                />
              </div>

              {/* Global parameters */}
              <div className="bg-gray-900/50 rounded-lg p-3 space-y-2 border border-gray-800">
                <div className="text-xs font-semibold text-gray-300 mb-2">Global Parameters</div>
                <div className="grid grid-cols-2 gap-x-4 gap-y-2">
                  <ParamInput
                    label="Tower Height"
                    value={tower.tower_height}
                    unit="m"
                    onChange={(v) => updateField('tower_height', v)}
                  />
                  <ParamInput
                    label="Base Height"
                    value={tower.tower_base_height}
                    unit="m"
                    onChange={(v) => updateField('tower_base_height', v)}
                  />
                  <ParamInput
                    label="FA Damp 1"
                    value={tower.tower_fa_damping_1}
                    unit="%"
                    onChange={(v) => updateField('tower_fa_damping_1', v)}
                  />
                  <ParamInput
                    label="FA Damp 2"
                    value={tower.tower_fa_damping_2}
                    unit="%"
                    onChange={(v) => updateField('tower_fa_damping_2', v)}
                  />
                  <ParamInput
                    label="SS Damp 1"
                    value={tower.tower_ss_damping_1}
                    unit="%"
                    onChange={(v) => updateField('tower_ss_damping_1', v)}
                  />
                  <ParamInput
                    label="SS Damp 2"
                    value={tower.tower_ss_damping_2}
                    unit="%"
                    onChange={(v) => updateField('tower_ss_damping_2', v)}
                  />
                </div>
              </div>

              {/* Station data grid */}
              <div className="text-xs font-semibold text-gray-300 mt-4 mb-1">
                Tower Stations
              </div>
              <DataGrid<TowerStation>
                columns={STATION_COLUMNS}
                rows={stations}
                onChange={handleStationsChange}
                onAddRow={handleAddStation}
                onDeleteRow={handleDeleteStation}
                maxHeight="calc(100vh - 450px)"
              />
            </div>
          </div>

          {/* CENTER PANEL: Visualization (~35%) */}
          <div className="w-[35%] flex flex-col border-r border-gray-800 overflow-hidden">
            <div className="px-3 py-2 bg-gray-900/50 border-b border-gray-800 shrink-0">
              <span className="text-xs font-semibold text-gray-300">Tower Profile</span>
            </div>
            <div className="flex-1 p-4 overflow-auto flex items-center justify-center">
              {stations.length > 0 ? (
                <TowerVisualization
                  stations={stations}
                  towerHeight={tower.tower_height}
                  towerBaseHeight={tower.tower_base_height}
                />
              ) : (
                <div className="text-center text-gray-500">
                  <Box size={48} className="mx-auto mb-3 opacity-30" />
                  <p className="text-sm">No stations defined</p>
                  <p className="text-xs mt-1">Add stations to see the tower profile</p>
                </div>
              )}
            </div>
          </div>

          {/* RIGHT PANEL: Charts + Preview (~25%) */}
          <div className="w-[25%] flex flex-col overflow-hidden">
            {/* Tabs */}
            <div className="flex border-b border-gray-800 shrink-0">
              <button
                onClick={() => setRightTab('plots')}
                className={clsx(
                  'flex-1 px-3 py-2 text-xs font-medium transition-colors',
                  rightTab === 'plots'
                    ? 'bg-gray-800 text-blue-400 border-b-2 border-blue-400'
                    : 'text-gray-400 hover:text-gray-200 hover:bg-gray-900'
                )}
              >
                Property Plots
              </button>
              <button
                onClick={() => setRightTab('preview')}
                className={clsx(
                  'flex-1 px-3 py-2 text-xs font-medium transition-colors',
                  rightTab === 'preview'
                    ? 'bg-gray-800 text-blue-400 border-b-2 border-blue-400'
                    : 'text-gray-400 hover:text-gray-200 hover:bg-gray-900'
                )}
              >
                File Preview
              </button>
            </div>

            <div className="flex-1 overflow-auto">
              {rightTab === 'plots' ? (
                <div className="p-2 space-y-2">
                  <PlotPanel
                    data={plotData.massPlot.data}
                    layout={plotData.massPlot.layout}
                    title={plotData.massPlot.title}
                    height={200}
                  />
                  <PlotPanel
                    data={plotData.stiffPlot.data}
                    layout={plotData.stiffPlot.layout}
                    title={plotData.stiffPlot.title}
                    height={200}
                  />
                  <PlotPanel
                    data={plotData.odPlot.data}
                    layout={plotData.odPlot.layout}
                    title={plotData.odPlot.title}
                    height={200}
                  />

                  {/* Mode shape coefficients */}
                  <div className="bg-gray-900/50 rounded-lg p-3 border border-gray-800 space-y-3">
                    <div className="text-xs font-semibold text-gray-300">
                      Mode Shape Coefficients
                    </div>
                    <ModeShapeInputs
                      label="FA Mode 1"
                      coeffs={tower.fa_mode_1_coeffs ?? [0, 0, 0, 0, 0]}
                      onChange={(c) => updateField('fa_mode_1_coeffs', c)}
                    />
                    <ModeShapeInputs
                      label="FA Mode 2"
                      coeffs={tower.fa_mode_2_coeffs ?? [0, 0, 0, 0, 0]}
                      onChange={(c) => updateField('fa_mode_2_coeffs', c)}
                    />
                    <ModeShapeInputs
                      label="SS Mode 1"
                      coeffs={tower.ss_mode_1_coeffs ?? [0, 0, 0, 0, 0]}
                      onChange={(c) => updateField('ss_mode_1_coeffs', c)}
                    />
                    <ModeShapeInputs
                      label="SS Mode 2"
                      coeffs={tower.ss_mode_2_coeffs ?? [0, 0, 0, 0, 0]}
                      onChange={(c) => updateField('ss_mode_2_coeffs', c)}
                    />
                  </div>
                </div>
              ) : (
                <div className="p-2 h-full">
                  <FilePreview
                    content={filePreview}
                    title="ElastoDyn Tower File"
                    isLoading={filePreviewLoading}
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
            <Box size={64} className="mx-auto mb-4 opacity-20" />
            <p className="text-lg font-medium">No tower selected</p>
            <p className="text-sm mt-2">
              Select a tower from the dropdown or create a new one
            </p>
            <button
              onClick={handleNewTower}
              className="mt-4 px-4 py-2 bg-blue-600 text-white rounded
                         hover:bg-blue-700 transition-colors text-sm"
            >
              Create New Tower
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
