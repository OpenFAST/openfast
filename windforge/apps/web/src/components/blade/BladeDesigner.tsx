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
  Wind,
} from 'lucide-react';
import clsx from 'clsx';
import DataGrid, { DataGridColumn } from '@/components/common/DataGrid';
import FilePreview from '@/components/common/FilePreview';
import PlotPanel from '@/components/common/PlotPanel';

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

interface BladeStructuralStation {
  frac: number;
  pitch_axis: number;
  struct_twist: number;
  mass_den: number;
  flap_stiff: number;
  edge_stiff: number;
}

interface BladeAeroStation {
  frac: number;
  chord: number;
  aero_twist: number;
  airfoil_id: string;
  aero_center: number;
}

interface Blade {
  id: string;
  project_id: string;
  name: string;
  version: number;
  blade_length: number;
  structural_stations: BladeStructuralStation[] | null;
  aero_stations: BladeAeroStation[] | null;
  flap_mode_1_coeffs: number[] | null;
  flap_mode_2_coeffs: number[] | null;
  edge_mode_1_coeffs: number[] | null;
  blade_flap_damping: number;
  blade_edge_damping: number;
  is_active: boolean;
  created_at: string;
}

const API_BASE = '/api/v1';

const STRUCTURAL_COLUMNS: DataGridColumn[] = [
  { name: 'Frac', key: 'frac', type: 'number', min: 0, max: 1, step: 0.01, unit: '-', precision: 4 },
  { name: 'Pitch Axis', key: 'pitch_axis', type: 'number', min: 0, max: 1, step: 0.01, unit: '-', precision: 4 },
  { name: 'Struct Twist', key: 'struct_twist', type: 'number', step: 0.1, unit: 'deg', precision: 3 },
  { name: 'Mass Den', key: 'mass_den', type: 'number', min: 0, step: 1, unit: 'kg/m', precision: 2 },
  { name: 'Flap Stiff', key: 'flap_stiff', type: 'number', min: 0, step: 1e6, unit: 'N-m^2', precision: 2 },
  { name: 'Edge Stiff', key: 'edge_stiff', type: 'number', min: 0, step: 1e6, unit: 'N-m^2', precision: 2 },
];

const AERO_COLUMNS: DataGridColumn[] = [
  { name: 'Frac', key: 'frac', type: 'number', min: 0, max: 1, step: 0.01, unit: '-', precision: 4 },
  { name: 'Chord', key: 'chord', type: 'number', min: 0, step: 0.01, unit: 'm', precision: 4 },
  { name: 'Aero Twist', key: 'aero_twist', type: 'number', step: 0.1, unit: 'deg', precision: 3 },
  { name: 'Airfoil', key: 'airfoil_id', type: 'text', unit: '' },
  { name: 'Aero Center', key: 'aero_center', type: 'number', min: 0, max: 1, step: 0.01, unit: '-', precision: 4 },
];

const DEFAULT_STRUCTURAL_STATION: BladeStructuralStation = {
  frac: 0,
  pitch_axis: 0.25,
  struct_twist: 0,
  mass_den: 0,
  flap_stiff: 0,
  edge_stiff: 0,
};

const DEFAULT_AERO_STATION: BladeAeroStation = {
  frac: 0,
  chord: 1.0,
  aero_twist: 0,
  airfoil_id: '',
  aero_center: 0.25,
};

// ---------------------------------------------------------------------------
// Helpers
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
      {unit && <span className="text-xs text-gray-500 w-12 shrink-0">{unit}</span>}
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
// Blade planform SVG
// ---------------------------------------------------------------------------

function BladePlanformView({
  aeroStations,
  bladeLength,
}: {
  aeroStations: BladeAeroStation[];
  bladeLength: number;
}) {
  const svgWidth = 300;
  const svgHeight = 500;
  const padTop = 30;
  const padBottom = 30;
  const padLeft = 50;
  const padRight = 30;
  const drawW = svgWidth - padLeft - padRight;
  const drawH = svgHeight - padTop - padBottom;

  const sorted = useMemo(
    () => [...aeroStations].sort((a, b) => a.frac - b.frac),
    [aeroStations]
  );

  const maxChord = useMemo(
    () => Math.max(...sorted.map((s) => s.chord), 1),
    [sorted]
  );

  // Blade root at bottom, tip at top
  // Chord envelope: LE on the left, TE on the right
  const lePoints: [number, number][] = [];
  const tePoints: [number, number][] = [];
  const twistLine: [number, number][] = [];

  const maxTwist = useMemo(
    () => Math.max(...sorted.map((s) => Math.abs(s.aero_twist)), 1),
    [sorted]
  );

  sorted.forEach((st) => {
    const y = padTop + drawH * (1 - st.frac);
    const chordPx = (st.chord / maxChord) * drawW * 0.7;
    const centerX = padLeft + drawW * 0.35;

    // Pitch axis offset
    const leX = centerX - chordPx * st.aero_center;
    const teX = centerX + chordPx * (1 - st.aero_center);

    lePoints.push([leX, y]);
    tePoints.push([teX, y]);

    // Twist overlay as secondary x position
    const twistX = padLeft + drawW * 0.85 + (st.aero_twist / maxTwist) * drawW * 0.12;
    twistLine.push([twistX, y]);
  });

  const chordPath =
    lePoints.length > 1
      ? [
          ...lePoints.map((p, i) => `${i === 0 ? 'M' : 'L'}${p[0]},${p[1]}`),
          ...tePoints.reverse().map((p) => `L${p[0]},${p[1]}`),
          'Z',
        ].join(' ')
      : '';

  const twistPath =
    twistLine.length > 1
      ? twistLine.map((p, i) => `${i === 0 ? 'M' : 'L'}${p[0]},${p[1]}`).join(' ')
      : '';

  // Span labels
  const spanLabels = useMemo(() => {
    const labels: { frac: number; span: number }[] = [];
    for (let f = 0; f <= 1.001; f += 0.2) {
      labels.push({ frac: Math.min(f, 1), span: Math.min(f, 1) * bladeLength });
    }
    return labels;
  }, [bladeLength]);

  return (
    <svg
      viewBox={`0 0 ${svgWidth} ${svgHeight}`}
      className="w-full h-full"
      style={{ maxHeight: '100%' }}
    >
      {/* Background grid */}
      <defs>
        <pattern id="bladeGrid" width="20" height="20" patternUnits="userSpaceOnUse">
          <path d="M 20 0 L 0 0 0 20" fill="none" stroke="rgba(75,85,99,0.15)" strokeWidth="0.5" />
        </pattern>
        <linearGradient id="bladeGradient" x1="0" y1="0" x2="1" y2="0">
          <stop offset="0%" stopColor="rgba(16,185,129,0.12)" />
          <stop offset="50%" stopColor="rgba(16,185,129,0.22)" />
          <stop offset="100%" stopColor="rgba(16,185,129,0.12)" />
        </linearGradient>
      </defs>
      <rect width={svgWidth} height={svgHeight} fill="url(#bladeGrid)" />

      {/* Root line */}
      <line
        x1={padLeft - 5}
        y1={padTop + drawH}
        x2={padLeft + drawW * 0.75}
        y2={padTop + drawH}
        stroke="#6b7280"
        strokeWidth={1.5}
      />
      <text
        x={padLeft + drawW * 0.35}
        y={padTop + drawH + 16}
        textAnchor="middle"
        fontSize={9}
        className="fill-gray-500"
      >
        Root
      </text>

      {/* Tip label */}
      <text
        x={padLeft + drawW * 0.35}
        y={padTop - 8}
        textAnchor="middle"
        fontSize={9}
        className="fill-gray-500"
      >
        Tip
      </text>

      {/* Chord envelope */}
      {chordPath && (
        <path d={chordPath} fill="url(#bladeGradient)" stroke="#10b981" strokeWidth={1.5} />
      )}

      {/* Twist line */}
      {twistPath && (
        <>
          {/* Twist axis */}
          <line
            x1={padLeft + drawW * 0.85}
            y1={padTop}
            x2={padLeft + drawW * 0.85}
            y2={padTop + drawH}
            stroke="rgba(245,158,11,0.2)"
            strokeWidth={0.5}
          />
          <path d={twistPath} fill="none" stroke="#f59e0b" strokeWidth={1.5} />
          <text
            x={padLeft + drawW * 0.85}
            y={padTop - 8}
            textAnchor="middle"
            fontSize={8}
            className="fill-amber-400"
          >
            Twist
          </text>
        </>
      )}

      {/* Station markers */}
      {sorted.map((st, i) => {
        const y = padTop + drawH * (1 - st.frac);
        const chordPx = (st.chord / maxChord) * drawW * 0.7;
        const centerX = padLeft + drawW * 0.35;
        const leX = centerX - chordPx * st.aero_center;
        return (
          <g key={i}>
            <circle cx={leX} cy={y} r={2.5} fill="#10b981" stroke="#064e3b" strokeWidth={0.8} />
            <circle
              cx={leX + chordPx}
              cy={y}
              r={2.5}
              fill="#10b981"
              stroke="#064e3b"
              strokeWidth={0.8}
            />
          </g>
        );
      })}

      {/* Span labels on left */}
      {spanLabels.map((sl, i) => {
        const y = padTop + drawH * (1 - sl.frac);
        return (
          <g key={i}>
            <line
              x1={8}
              y1={y}
              x2={padLeft - 8}
              y2={y}
              stroke="rgba(156,163,175,0.2)"
              strokeWidth={0.5}
            />
            <text x={5} y={y + 3} fontSize={8} className="fill-gray-500" textAnchor="start">
              {sl.span.toFixed(1)}
            </text>
          </g>
        );
      })}
      <text x={8} y={padTop - 12} fontSize={9} className="fill-gray-400 font-semibold">
        Span (m)
      </text>
    </svg>
  );
}

// ---------------------------------------------------------------------------
// Main BladeDesigner component
// ---------------------------------------------------------------------------

export default function BladeDesigner() {
  const { projectId } = useParams<{ projectId: string }>();

  // State
  const [blades, setBlades] = useState<Blade[]>([]);
  const [selectedBladeId, setSelectedBladeId] = useState<string | null>(null);
  const [blade, setBlade] = useState<Blade | null>(null);
  const [originalBlade, setOriginalBlade] = useState<Blade | null>(null);
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [leftTab, setLeftTab] = useState<'structural' | 'aero'>('structural');
  const [rightTab, setRightTab] = useState<'distributions' | 'ed_preview' | 'ad_preview'>(
    'distributions'
  );
  const [edPreview, setEdPreview] = useState<string | null>(null);
  const [adPreview, setAdPreview] = useState<string | null>(null);
  const [previewLoading, setPreviewLoading] = useState(false);
  const [dropdownOpen, setDropdownOpen] = useState(false);

  const structuralStations = blade?.structural_stations ?? [];
  const aeroStations = blade?.aero_stations ?? [];

  const isDirty = useMemo(
    () => JSON.stringify(blade) !== JSON.stringify(originalBlade),
    [blade, originalBlade]
  );

  // ---------------------------------------------------------------------------
  // Fetching
  // ---------------------------------------------------------------------------

  const fetchBlades = useCallback(async () => {
    if (!projectId) return;
    try {
      setLoading(true);
      const resp = await axios.get<Blade[]>(`${API_BASE}/projects/${projectId}/blades`);
      setBlades(resp.data);
      if (resp.data.length > 0 && !selectedBladeId) {
        setSelectedBladeId(resp.data[0].id);
      }
    } catch (err: any) {
      setError(err.message || 'Failed to load blades');
    } finally {
      setLoading(false);
    }
  }, [projectId, selectedBladeId]);

  const fetchBlade = useCallback(
    async (bladeId: string) => {
      if (!projectId) return;
      try {
        setLoading(true);
        const resp = await axios.get<Blade>(
          `${API_BASE}/projects/${projectId}/blades/${bladeId}`
        );
        setBlade(resp.data);
        setOriginalBlade(JSON.parse(JSON.stringify(resp.data)));
        setError(null);
      } catch (err: any) {
        setError(err.message || 'Failed to load blade');
      } finally {
        setLoading(false);
      }
    },
    [projectId]
  );

  const fetchEdPreview = useCallback(async () => {
    if (!projectId || !selectedBladeId) return;
    try {
      setPreviewLoading(true);
      const resp = await axios.get<string>(
        `${API_BASE}/projects/${projectId}/blades/${selectedBladeId}/preview-ed`
      );
      setEdPreview(resp.data);
    } catch {
      setEdPreview(null);
    } finally {
      setPreviewLoading(false);
    }
  }, [projectId, selectedBladeId]);

  const fetchAdPreview = useCallback(async () => {
    if (!projectId || !selectedBladeId) return;
    try {
      setPreviewLoading(true);
      const resp = await axios.get<string>(
        `${API_BASE}/projects/${projectId}/blades/${selectedBladeId}/preview-ad`
      );
      setAdPreview(resp.data);
    } catch {
      setAdPreview(null);
    } finally {
      setPreviewLoading(false);
    }
  }, [projectId, selectedBladeId]);

  useEffect(() => {
    fetchBlades();
  }, [fetchBlades]);

  useEffect(() => {
    if (selectedBladeId) fetchBlade(selectedBladeId);
  }, [selectedBladeId, fetchBlade]);

  useEffect(() => {
    if (rightTab === 'ed_preview') fetchEdPreview();
    if (rightTab === 'ad_preview') fetchAdPreview();
  }, [rightTab, fetchEdPreview, fetchAdPreview]);

  // ---------------------------------------------------------------------------
  // Actions
  // ---------------------------------------------------------------------------

  const handleSave = useCallback(async () => {
    if (!projectId || !blade) return;
    try {
      setSaving(true);
      if (blade.id) {
        await axios.put(
          `${API_BASE}/projects/${projectId}/blades/${blade.id}`,
          blade
        );
      } else {
        const resp = await axios.post<Blade>(
          `${API_BASE}/projects/${projectId}/blades`,
          blade
        );
        setSelectedBladeId(resp.data.id);
      }
      setOriginalBlade(JSON.parse(JSON.stringify(blade)));
      await fetchBlades();
      setError(null);
    } catch (err: any) {
      setError(err.message || 'Failed to save blade');
    } finally {
      setSaving(false);
    }
  }, [projectId, blade, fetchBlades]);

  const handleRevert = useCallback(() => {
    if (originalBlade) setBlade(JSON.parse(JSON.stringify(originalBlade)));
  }, [originalBlade]);

  const handleNewBlade = useCallback(() => {
    const newBlade: Blade = {
      id: '',
      project_id: projectId || '',
      name: 'New Blade',
      version: 1,
      blade_length: 61.5,
      structural_stations: [
        { frac: 0.0, pitch_axis: 0.25, struct_twist: 13.308, mass_den: 678.935, flap_stiff: 1.81e10, edge_stiff: 1.81e10 },
        { frac: 0.25, pitch_axis: 0.25, struct_twist: 9.0, mass_den: 450.0, flap_stiff: 9.0e9, edge_stiff: 1.0e10 },
        { frac: 0.5, pitch_axis: 0.25, struct_twist: 4.0, mass_den: 350.0, flap_stiff: 4.0e9, edge_stiff: 5.5e9 },
        { frac: 0.75, pitch_axis: 0.375, struct_twist: 1.0, mass_den: 200.0, flap_stiff: 1.0e9, edge_stiff: 1.5e9 },
        { frac: 1.0, pitch_axis: 0.5, struct_twist: 0.0, mass_den: 10.0, flap_stiff: 1.0e7, edge_stiff: 1.0e7 },
      ],
      aero_stations: [
        { frac: 0.0, chord: 3.542, aero_twist: 13.308, airfoil_id: 'cylinder', aero_center: 0.25 },
        { frac: 0.15, chord: 4.557, aero_twist: 11.48, airfoil_id: 'DU40', aero_center: 0.25 },
        { frac: 0.35, chord: 4.167, aero_twist: 6.544, airfoil_id: 'DU30', aero_center: 0.25 },
        { frac: 0.55, chord: 3.502, aero_twist: 2.87, airfoil_id: 'DU25', aero_center: 0.25 },
        { frac: 0.75, chord: 2.764, aero_twist: 0.724, airfoil_id: 'DU21', aero_center: 0.25 },
        { frac: 1.0, chord: 1.419, aero_twist: -1.0, airfoil_id: 'NACA64', aero_center: 0.25 },
      ],
      flap_mode_1_coeffs: [0.0622, 1.7254, -3.2452, 4.7131, -2.2555],
      flap_mode_2_coeffs: [-0.5809, 1.2067, -15.5349, 29.7347, -13.8245],
      edge_mode_1_coeffs: [0.5994, -0.2658, -4.0018, 7.6563, -2.9888],
      blade_flap_damping: 2.5,
      blade_edge_damping: 2.5,
      is_active: true,
      created_at: new Date().toISOString(),
    };
    setBlade(newBlade);
    setOriginalBlade(null);
    setSelectedBladeId(null);
  }, [projectId]);

  const handleAddStructuralStation = useCallback(() => {
    if (!blade) return;
    const current = blade.structural_stations ?? [];
    const lastFrac = current.length > 0 ? current[current.length - 1].frac : 0;
    setBlade({
      ...blade,
      structural_stations: [
        ...current,
        { ...DEFAULT_STRUCTURAL_STATION, frac: Math.min(lastFrac + 0.05, 1.0) },
      ],
    });
  }, [blade]);

  const handleDeleteStructuralStation = useCallback(
    (index: number) => {
      if (!blade || !blade.structural_stations) return;
      setBlade({
        ...blade,
        structural_stations: blade.structural_stations.filter((_, i) => i !== index),
      });
    },
    [blade]
  );

  const handleAddAeroStation = useCallback(() => {
    if (!blade) return;
    const current = blade.aero_stations ?? [];
    const lastFrac = current.length > 0 ? current[current.length - 1].frac : 0;
    setBlade({
      ...blade,
      aero_stations: [
        ...current,
        { ...DEFAULT_AERO_STATION, frac: Math.min(lastFrac + 0.05, 1.0) },
      ],
    });
  }, [blade]);

  const handleDeleteAeroStation = useCallback(
    (index: number) => {
      if (!blade || !blade.aero_stations) return;
      setBlade({
        ...blade,
        aero_stations: blade.aero_stations.filter((_, i) => i !== index),
      });
    },
    [blade]
  );

  const updateField = useCallback(
    <K extends keyof Blade>(key: K, value: Blade[K]) => {
      if (!blade) return;
      setBlade({ ...blade, [key]: value });
    },
    [blade]
  );

  // ---------------------------------------------------------------------------
  // Plot data
  // ---------------------------------------------------------------------------

  const plotData = useMemo(() => {
    const sortedStructural = [...structuralStations].sort((a, b) => a.frac - b.frac);
    const sortedAero = [...aeroStations].sort((a, b) => a.frac - b.frac);
    const bladeLen = blade?.blade_length ?? 1;

    const chordPlot = {
      data: [
        {
          x: sortedAero.map((s) => s.frac * bladeLen),
          y: sortedAero.map((s) => s.chord),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'Chord',
          line: { color: '#10b981', width: 2 },
          marker: { size: 4 },
        },
      ],
      layout: {
        xaxis: { title: { text: 'Span (m)' } },
        yaxis: { title: { text: 'm' } },
      },
      title: 'Chord Distribution',
    };

    const twistPlot = {
      data: [
        {
          x: sortedAero.map((s) => s.frac * bladeLen),
          y: sortedAero.map((s) => s.aero_twist),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'Aero Twist',
          line: { color: '#f59e0b', width: 2 },
          marker: { size: 4 },
        },
        {
          x: sortedStructural.map((s) => s.frac * bladeLen),
          y: sortedStructural.map((s) => s.struct_twist),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'Structural Twist',
          line: { color: '#ef4444', width: 2, dash: 'dash' as const },
          marker: { size: 4 },
        },
      ],
      layout: {
        xaxis: { title: { text: 'Span (m)' } },
        yaxis: { title: { text: 'deg' } },
      },
      title: 'Twist Distribution',
    };

    const massPlot = {
      data: [
        {
          x: sortedStructural.map((s) => s.frac * bladeLen),
          y: sortedStructural.map((s) => s.mass_den),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'Mass Density',
          line: { color: '#3b82f6', width: 2 },
          marker: { size: 4 },
        },
      ],
      layout: {
        xaxis: { title: { text: 'Span (m)' } },
        yaxis: { title: { text: 'kg/m' } },
      },
      title: 'Mass Distribution',
    };

    const stiffPlot = {
      data: [
        {
          x: sortedStructural.map((s) => s.frac * bladeLen),
          y: sortedStructural.map((s) => s.flap_stiff),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'Flap Stiffness',
          line: { color: '#8b5cf6', width: 2 },
          marker: { size: 4 },
        },
        {
          x: sortedStructural.map((s) => s.frac * bladeLen),
          y: sortedStructural.map((s) => s.edge_stiff),
          type: 'scatter' as const,
          mode: 'lines+markers' as const,
          name: 'Edge Stiffness',
          line: { color: '#ec4899', width: 2, dash: 'dash' as const },
          marker: { size: 4 },
        },
      ],
      layout: {
        xaxis: { title: { text: 'Span (m)' } },
        yaxis: { title: { text: 'N-m^2' } },
      },
      title: 'Stiffness Distribution',
    };

    return { chordPlot, twistPlot, massPlot, stiffPlot };
  }, [structuralStations, aeroStations, blade?.blade_length]);

  // ---------------------------------------------------------------------------
  // Render
  // ---------------------------------------------------------------------------

  if (loading && !blade) {
    return (
      <div className="flex items-center justify-center h-full bg-gray-950 text-gray-400">
        <Loader2 className="animate-spin mr-2" size={20} />
        Loading blade data...
      </div>
    );
  }

  return (
    <div className="flex flex-col h-full bg-gray-950 text-gray-200">
      {/* Top bar */}
      <div className="flex items-center gap-3 px-4 py-2 bg-gray-900 border-b border-gray-800 shrink-0">
        <Wind size={16} className="text-green-400" />
        <span className="text-sm font-semibold text-gray-200">Blade Designer</span>

        {/* Blade selector */}
        <div className="relative ml-4">
          <button
            onClick={() => setDropdownOpen(!dropdownOpen)}
            className="flex items-center gap-2 px-3 py-1.5 bg-gray-800 border border-gray-700
                       rounded text-xs hover:bg-gray-750 transition-colors min-w-[200px]"
          >
            <span className="truncate">{blade?.name || 'Select Blade'}</span>
            <ChevronDown size={14} className="ml-auto text-gray-500" />
          </button>
          {dropdownOpen && (
            <>
              <div className="fixed inset-0 z-10" onClick={() => setDropdownOpen(false)} />
              <div className="absolute top-full left-0 mt-1 w-64 bg-gray-800 border border-gray-700
                              rounded-lg shadow-xl z-20 max-h-60 overflow-auto">
                {blades.map((b) => (
                  <button
                    key={b.id}
                    onClick={() => {
                      setSelectedBladeId(b.id);
                      setDropdownOpen(false);
                    }}
                    className={clsx(
                      'w-full text-left px-3 py-2 text-xs hover:bg-gray-700 transition-colors',
                      b.id === selectedBladeId && 'bg-green-900/40 text-green-300'
                    )}
                  >
                    <div className="font-medium">{b.name}</div>
                    <div className="text-gray-500">v{b.version} &middot; {b.blade_length}m</div>
                  </button>
                ))}
                {blades.length === 0 && (
                  <div className="px-3 py-4 text-xs text-gray-500 text-center">No blades yet</div>
                )}
              </div>
            </>
          )}
        </div>

        <button
          onClick={handleNewBlade}
          className="flex items-center gap-1 px-3 py-1.5 text-xs bg-gray-800
                     border border-gray-700 rounded hover:bg-gray-750 transition-colors"
        >
          <Plus size={14} />
          New Blade
        </button>

        {blade && (
          <span className="ml-2 px-2 py-0.5 bg-green-900/40 text-green-300 text-xs rounded-full font-medium">
            v{blade.version}
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
            disabled={saving || !blade}
            className="flex items-center gap-1 px-4 py-1.5 text-xs bg-green-600
                       text-white rounded hover:bg-green-700
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

      {/* Main content */}
      {blade ? (
        <div className="flex-1 flex overflow-hidden">
          {/* LEFT PANEL: Data (~40%) */}
          <div className="w-[40%] flex flex-col border-r border-gray-800 overflow-hidden">
            <div className="p-3 space-y-3 overflow-auto flex-1">
              {/* Blade name */}
              <div className="flex items-center gap-2">
                <label className="text-xs text-gray-400 w-20 shrink-0">Name</label>
                <input
                  type="text"
                  value={blade.name}
                  onChange={(e) => updateField('name', e.target.value)}
                  className="flex-1 px-2 py-1.5 bg-gray-800 border border-gray-700 rounded
                             text-sm text-gray-200 focus:border-green-500
                             focus:outline-none focus:ring-1 focus:ring-green-500/50"
                />
              </div>

              {/* Global params */}
              <div className="bg-gray-900/50 rounded-lg p-3 space-y-2 border border-gray-800">
                <div className="text-xs font-semibold text-gray-300 mb-2">Global Parameters</div>
                <ParamInput
                  label="Blade Length"
                  value={blade.blade_length}
                  unit="m"
                  onChange={(v) => updateField('blade_length', v)}
                />
                <ParamInput
                  label="Flap Damping"
                  value={blade.blade_flap_damping}
                  unit="%"
                  onChange={(v) => updateField('blade_flap_damping', v)}
                />
                <ParamInput
                  label="Edge Damping"
                  value={blade.blade_edge_damping}
                  unit="%"
                  onChange={(v) => updateField('blade_edge_damping', v)}
                />
              </div>

              {/* Tabs: Structural / Aero */}
              <div className="flex border-b border-gray-800">
                <button
                  onClick={() => setLeftTab('structural')}
                  className={clsx(
                    'flex-1 px-3 py-2 text-xs font-medium transition-colors',
                    leftTab === 'structural'
                      ? 'text-green-400 border-b-2 border-green-400'
                      : 'text-gray-400 hover:text-gray-200'
                  )}
                >
                  Structural ({structuralStations.length})
                </button>
                <button
                  onClick={() => setLeftTab('aero')}
                  className={clsx(
                    'flex-1 px-3 py-2 text-xs font-medium transition-colors',
                    leftTab === 'aero'
                      ? 'text-amber-400 border-b-2 border-amber-400'
                      : 'text-gray-400 hover:text-gray-200'
                  )}
                >
                  Aerodynamic ({aeroStations.length})
                </button>
              </div>

              {leftTab === 'structural' ? (
                <DataGrid<BladeStructuralStation>
                  columns={STRUCTURAL_COLUMNS}
                  rows={structuralStations}
                  onChange={(rows) => updateField('structural_stations', rows)}
                  onAddRow={handleAddStructuralStation}
                  onDeleteRow={handleDeleteStructuralStation}
                  maxHeight="calc(100vh - 480px)"
                />
              ) : (
                <DataGrid<BladeAeroStation>
                  columns={AERO_COLUMNS}
                  rows={aeroStations}
                  onChange={(rows) => updateField('aero_stations', rows)}
                  onAddRow={handleAddAeroStation}
                  onDeleteRow={handleDeleteAeroStation}
                  maxHeight="calc(100vh - 480px)"
                />
              )}
            </div>
          </div>

          {/* CENTER PANEL: Visualization (~35%) */}
          <div className="w-[35%] flex flex-col border-r border-gray-800 overflow-hidden">
            <div className="px-3 py-2 bg-gray-900/50 border-b border-gray-800 shrink-0">
              <span className="text-xs font-semibold text-gray-300">Blade Planform</span>
            </div>
            <div className="flex-1 p-4 overflow-auto flex items-center justify-center">
              {aeroStations.length > 0 ? (
                <BladePlanformView
                  aeroStations={aeroStations}
                  bladeLength={blade.blade_length}
                />
              ) : (
                <div className="text-center text-gray-500">
                  <Wind size={48} className="mx-auto mb-3 opacity-30" />
                  <p className="text-sm">No aero stations defined</p>
                  <p className="text-xs mt-1">Add stations to see the blade planform</p>
                </div>
              )}
            </div>
          </div>

          {/* RIGHT PANEL: Charts + Preview (~25%) */}
          <div className="w-[25%] flex flex-col overflow-hidden">
            {/* Tabs */}
            <div className="flex border-b border-gray-800 shrink-0">
              {(['distributions', 'ed_preview', 'ad_preview'] as const).map((tab) => (
                <button
                  key={tab}
                  onClick={() => setRightTab(tab)}
                  className={clsx(
                    'flex-1 px-2 py-2 text-xs font-medium transition-colors whitespace-nowrap',
                    rightTab === tab
                      ? 'bg-gray-800 text-green-400 border-b-2 border-green-400'
                      : 'text-gray-400 hover:text-gray-200 hover:bg-gray-900'
                  )}
                >
                  {tab === 'distributions'
                    ? 'Distributions'
                    : tab === 'ed_preview'
                    ? 'ED Preview'
                    : 'AD Preview'}
                </button>
              ))}
            </div>

            <div className="flex-1 overflow-auto">
              {rightTab === 'distributions' ? (
                <div className="p-2 space-y-2">
                  <PlotPanel
                    data={plotData.chordPlot.data}
                    layout={plotData.chordPlot.layout}
                    title={plotData.chordPlot.title}
                    height={180}
                  />
                  <PlotPanel
                    data={plotData.twistPlot.data}
                    layout={plotData.twistPlot.layout}
                    title={plotData.twistPlot.title}
                    height={180}
                  />
                  <PlotPanel
                    data={plotData.massPlot.data}
                    layout={plotData.massPlot.layout}
                    title={plotData.massPlot.title}
                    height={180}
                  />
                  <PlotPanel
                    data={plotData.stiffPlot.data}
                    layout={plotData.stiffPlot.layout}
                    title={plotData.stiffPlot.title}
                    height={180}
                  />

                  {/* Mode shape coefficients */}
                  <div className="bg-gray-900/50 rounded-lg p-3 border border-gray-800 space-y-3">
                    <div className="text-xs font-semibold text-gray-300">Mode Shape Coefficients</div>
                    <ModeShapeInputs
                      label="Flap Mode 1"
                      coeffs={blade.flap_mode_1_coeffs ?? [0, 0, 0, 0, 0]}
                      onChange={(c) => updateField('flap_mode_1_coeffs', c)}
                    />
                    <ModeShapeInputs
                      label="Flap Mode 2"
                      coeffs={blade.flap_mode_2_coeffs ?? [0, 0, 0, 0, 0]}
                      onChange={(c) => updateField('flap_mode_2_coeffs', c)}
                    />
                    <ModeShapeInputs
                      label="Edge Mode 1"
                      coeffs={blade.edge_mode_1_coeffs ?? [0, 0, 0, 0, 0]}
                      onChange={(c) => updateField('edge_mode_1_coeffs', c)}
                    />
                  </div>
                </div>
              ) : rightTab === 'ed_preview' ? (
                <div className="p-2 h-full">
                  <FilePreview
                    content={edPreview}
                    title="ElastoDyn Blade File"
                    isLoading={previewLoading}
                    className="h-full"
                  />
                </div>
              ) : (
                <div className="p-2 h-full">
                  <FilePreview
                    content={adPreview}
                    title="AeroDyn Blade File"
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
            <Wind size={64} className="mx-auto mb-4 opacity-20" />
            <p className="text-lg font-medium">No blade selected</p>
            <p className="text-sm mt-2">Select a blade from the dropdown or create a new one</p>
            <button
              onClick={handleNewBlade}
              className="mt-4 px-4 py-2 bg-green-600 text-white rounded
                         hover:bg-green-700 transition-colors text-sm"
            >
              Create New Blade
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
