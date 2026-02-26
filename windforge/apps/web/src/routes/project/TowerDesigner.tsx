import { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { towersApi } from '@/api/client';
import type { Tower } from '@/types';
import { Building2, Plus, Loader2, Trash2, ChevronDown } from 'lucide-react';
import clsx from 'clsx';

export default function TowerDesigner() {
  const { projectId } = useParams<{ projectId: string }>();
  const [towers, setTowers] = useState<Tower[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedTower, setSelectedTower] = useState<Tower | null>(null);

  useEffect(() => {
    if (!projectId) return;
    setIsLoading(true);
    towersApi
      .list(projectId)
      .then((data) => {
        setTowers(data);
        if (data.length > 0) setSelectedTower(data[0]);
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
          <h2 className="text-xl font-bold text-slate-800">Tower Designer</h2>
          <p className="text-sm text-slate-500">
            Define tower geometry, materials, and structural stations
          </p>
        </div>
        <button className="btn-primary">
          <Plus className="h-4 w-4" />
          New Tower
        </button>
      </div>

      {towers.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Building2 className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">
            No towers defined
          </h3>
          <p className="mt-1 text-sm text-slate-500">
            Create your first tower configuration to get started.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Tower list */}
          <div className="space-y-2">
            {towers.map((tower) => (
              <button
                key={tower.id}
                onClick={() => setSelectedTower(tower)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedTower?.id === tower.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <div className="flex items-center justify-between">
                  <span className="font-medium text-slate-800">
                    {tower.name}
                  </span>
                  <ChevronDown
                    className={clsx(
                      'h-4 w-4 text-slate-400 transition-transform',
                      selectedTower?.id === tower.id && 'rotate-180',
                    )}
                  />
                </div>
                <div className="mt-1 text-xs text-slate-500">
                  {tower.total_height_m}m tall &middot; {tower.num_stations}{' '}
                  stations
                </div>
              </button>
            ))}
          </div>

          {/* Tower detail */}
          {selectedTower && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-4">
                <h3 className="text-lg font-semibold text-slate-800">
                  {selectedTower.name}
                </h3>
                <button className="text-slate-400 hover:text-danger-500 transition-colors">
                  <Trash2 className="h-4 w-4" />
                </button>
              </div>

              <div className="grid grid-cols-3 gap-4 mb-6">
                <div className="rounded-lg bg-slate-50 p-3">
                  <p className="text-xs text-slate-500">Total Height</p>
                  <p className="text-lg font-semibold text-slate-800">
                    {selectedTower.total_height_m} m
                  </p>
                </div>
                <div className="rounded-lg bg-slate-50 p-3">
                  <p className="text-xs text-slate-500">Base Diameter</p>
                  <p className="text-lg font-semibold text-slate-800">
                    {selectedTower.base_diameter_m} m
                  </p>
                </div>
                <div className="rounded-lg bg-slate-50 p-3">
                  <p className="text-xs text-slate-500">Top Diameter</p>
                  <p className="text-lg font-semibold text-slate-800">
                    {selectedTower.top_diameter_m} m
                  </p>
                </div>
              </div>

              {/* Stations table */}
              {selectedTower.stations.length > 0 && (
                <div>
                  <h4 className="text-sm font-medium text-slate-700 mb-2">
                    Structural Stations
                  </h4>
                  <div className="overflow-x-auto rounded-lg border border-slate-200">
                    <table className="min-w-full text-sm">
                      <thead className="bg-slate-50">
                        <tr>
                          <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                            Height (m)
                          </th>
                          <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                            OD (m)
                          </th>
                          <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                            Thickness (m)
                          </th>
                          <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                            Density (kg/m3)
                          </th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-slate-100">
                        {selectedTower.stations.map((s, i) => (
                          <tr key={i} className="hover:bg-slate-50">
                            <td className="px-3 py-2 text-slate-700">
                              {s.height_m}
                            </td>
                            <td className="px-3 py-2 text-slate-700">
                              {s.outer_diameter_m}
                            </td>
                            <td className="px-3 py-2 text-slate-700">
                              {s.wall_thickness_m}
                            </td>
                            <td className="px-3 py-2 text-slate-700">
                              {s.density_kg_m3}
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
