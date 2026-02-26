import { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { bladesApi } from '@/api/client';
import type { Blade } from '@/types';
import { Fan, Plus, Loader2, Trash2, ChevronDown } from 'lucide-react';
import clsx from 'clsx';

export default function BladeDesigner() {
  const { projectId } = useParams<{ projectId: string }>();
  const [blades, setBlades] = useState<Blade[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedBlade, setSelectedBlade] = useState<Blade | null>(null);
  const [activeTab, setActiveTab] = useState<'structural' | 'aero'>(
    'structural',
  );

  useEffect(() => {
    if (!projectId) return;
    setIsLoading(true);
    bladesApi
      .list(projectId)
      .then((data) => {
        setBlades(data);
        if (data.length > 0) setSelectedBlade(data[0]);
      })
      .catch(() => {})
      .finally(() => setIsLoading(false));
  }, [projectId]);

  if (isLoading) {
    return (
      <div className="flex items-center justify-center py-24">
        <Loader2 className="h-8 w-8 animate-spin text-accent-500" />
      </div>
    );
  }

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-xl font-bold text-slate-800">Blade Designer</h2>
          <p className="text-sm text-slate-500">
            Define blade structural and aerodynamic properties
          </p>
        </div>
        <button className="btn-primary">
          <Plus className="h-4 w-4" />
          New Blade
        </button>
      </div>

      {blades.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Fan className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">
            No blades defined
          </h3>
          <p className="mt-1 text-sm text-slate-500">
            Create your first blade design to define structural and
            aerodynamic properties.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Blade list */}
          <div className="space-y-2">
            {blades.map((blade) => (
              <button
                key={blade.id}
                onClick={() => setSelectedBlade(blade)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedBlade?.id === blade.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <div className="flex items-center justify-between">
                  <span className="font-medium text-slate-800">
                    {blade.name}
                  </span>
                  <ChevronDown
                    className={clsx(
                      'h-4 w-4 text-slate-400 transition-transform',
                      selectedBlade?.id === blade.id && 'rotate-180',
                    )}
                  />
                </div>
                <div className="mt-1 text-xs text-slate-500">
                  {blade.length_m}m long &middot;{' '}
                  {blade.num_aero_stations} aero stations
                </div>
              </button>
            ))}
          </div>

          {/* Blade detail */}
          {selectedBlade && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-4">
                <h3 className="text-lg font-semibold text-slate-800">
                  {selectedBlade.name}
                </h3>
                <button className="text-slate-400 hover:text-danger-500 transition-colors">
                  <Trash2 className="h-4 w-4" />
                </button>
              </div>

              <div className="grid grid-cols-3 gap-4 mb-6">
                <div className="rounded-lg bg-slate-50 p-3">
                  <p className="text-xs text-slate-500">Blade Length</p>
                  <p className="text-lg font-semibold text-slate-800">
                    {selectedBlade.length_m} m
                  </p>
                </div>
                <div className="rounded-lg bg-slate-50 p-3">
                  <p className="text-xs text-slate-500">
                    Structural Stations
                  </p>
                  <p className="text-lg font-semibold text-slate-800">
                    {selectedBlade.num_structural_stations}
                  </p>
                </div>
                <div className="rounded-lg bg-slate-50 p-3">
                  <p className="text-xs text-slate-500">Aero Stations</p>
                  <p className="text-lg font-semibold text-slate-800">
                    {selectedBlade.num_aero_stations}
                  </p>
                </div>
              </div>

              {/* Tab toggle */}
              <div className="flex gap-1 rounded-lg bg-slate-100 p-1 mb-4">
                <button
                  onClick={() => setActiveTab('structural')}
                  className={clsx(
                    'flex-1 rounded-md px-3 py-1.5 text-sm font-medium transition-all',
                    activeTab === 'structural'
                      ? 'bg-white text-slate-800 shadow-sm'
                      : 'text-slate-500 hover:text-slate-700',
                  )}
                >
                  Structural
                </button>
                <button
                  onClick={() => setActiveTab('aero')}
                  className={clsx(
                    'flex-1 rounded-md px-3 py-1.5 text-sm font-medium transition-all',
                    activeTab === 'aero'
                      ? 'bg-white text-slate-800 shadow-sm'
                      : 'text-slate-500 hover:text-slate-700',
                  )}
                >
                  Aerodynamic
                </button>
              </div>

              {/* Tables */}
              <div className="overflow-x-auto rounded-lg border border-slate-200">
                {activeTab === 'structural' ? (
                  <table className="min-w-full text-sm">
                    <thead className="bg-slate-50">
                      <tr>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Span
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Mass (kg/m)
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Flap EI
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Edge EI
                        </th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-100">
                      {selectedBlade.structural_stations.map((s, i) => (
                        <tr key={i} className="hover:bg-slate-50">
                          <td className="px-3 py-2 text-slate-700">
                            {s.span_fraction.toFixed(3)}
                          </td>
                          <td className="px-3 py-2 text-slate-700">
                            {s.mass_density_kg_m.toFixed(1)}
                          </td>
                          <td className="px-3 py-2 text-slate-700">
                            {s.flap_stiffness_nm2.toExponential(2)}
                          </td>
                          <td className="px-3 py-2 text-slate-700">
                            {s.edge_stiffness_nm2.toExponential(2)}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                ) : (
                  <table className="min-w-full text-sm">
                    <thead className="bg-slate-50">
                      <tr>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Span
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Chord (m)
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Twist (deg)
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Airfoil
                        </th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-100">
                      {selectedBlade.aero_stations.map((s, i) => (
                        <tr key={i} className="hover:bg-slate-50">
                          <td className="px-3 py-2 text-slate-700">
                            {s.span_fraction.toFixed(3)}
                          </td>
                          <td className="px-3 py-2 text-slate-700">
                            {s.chord_m.toFixed(3)}
                          </td>
                          <td className="px-3 py-2 text-slate-700">
                            {s.twist_deg.toFixed(2)}
                          </td>
                          <td className="px-3 py-2 text-slate-700">
                            {s.airfoil_name}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                )}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
