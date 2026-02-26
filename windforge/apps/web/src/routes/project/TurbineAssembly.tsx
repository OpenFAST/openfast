import { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { turbineModelsApi, towersApi, bladesApi, controllersApi } from '@/api/client';
import type { TurbineModel, Tower, Blade, Controller } from '@/types';
import { Boxes, Plus, Loader2, Trash2 } from 'lucide-react';
import clsx from 'clsx';

export default function TurbineAssembly() {
  const { projectId } = useParams<{ projectId: string }>();
  const [models, setModels] = useState<TurbineModel[]>([]);
  const [towers, setTowers] = useState<Tower[]>([]);
  const [blades, setBlades] = useState<Blade[]>([]);
  const [controllers, setControllers] = useState<Controller[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedModel, setSelectedModel] = useState<TurbineModel | null>(null);

  useEffect(() => {
    if (!projectId) return;
    setIsLoading(true);
    Promise.all([
      turbineModelsApi.list(projectId),
      towersApi.list(projectId),
      bladesApi.list(projectId),
      controllersApi.list(projectId),
    ])
      .then(([modelsData, towersData, bladesData, controllersData]) => {
        setModels(modelsData);
        setTowers(towersData);
        setBlades(bladesData);
        setControllers(controllersData);
        if (modelsData.length > 0) setSelectedModel(modelsData[0]);
      })
      .catch(() => {})
      .finally(() => setIsLoading(false));
  }, [projectId]);

  const getTowerName = (id: string) =>
    towers.find((t) => t.id === id)?.name ?? 'Unknown';
  const getBladeName = (id: string) =>
    blades.find((b) => b.id === id)?.name ?? 'Unknown';
  const getControllerName = (id: string) =>
    controllers.find((c) => c.id === id)?.name ?? 'Unknown';

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
          <h2 className="text-xl font-bold text-slate-800">
            Turbine Assembly
          </h2>
          <p className="text-sm text-slate-500">
            Combine tower, blade, and controller into a complete turbine model
          </p>
        </div>
        <button className="btn-primary">
          <Plus className="h-4 w-4" />
          New Assembly
        </button>
      </div>

      {models.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Boxes className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">
            No turbine assemblies
          </h3>
          <p className="mt-1 text-sm text-slate-500 max-w-sm text-center">
            Define tower, blade, and controller first, then assemble them into a
            complete turbine model.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Model list */}
          <div className="space-y-2">
            {models.map((model) => (
              <button
                key={model.id}
                onClick={() => setSelectedModel(model)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedModel?.id === model.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <span className="font-medium text-slate-800">
                  {model.name}
                </span>
                <div className="mt-1 text-xs text-slate-500">
                  {model.num_blades} blades &middot; GR{' '}
                  {model.gearbox_ratio}
                </div>
              </button>
            ))}
          </div>

          {/* Model detail */}
          {selectedModel && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-6">
                <h3 className="text-lg font-semibold text-slate-800">
                  {selectedModel.name}
                </h3>
                <button className="text-slate-400 hover:text-danger-500 transition-colors">
                  <Trash2 className="h-4 w-4" />
                </button>
              </div>

              {/* Component references */}
              <div className="grid grid-cols-3 gap-4 mb-6">
                <div className="rounded-lg border border-slate-200 p-4">
                  <p className="text-xs text-slate-500 mb-1">Tower</p>
                  <p className="text-sm font-semibold text-slate-800">
                    {getTowerName(selectedModel.tower_id)}
                  </p>
                </div>
                <div className="rounded-lg border border-slate-200 p-4">
                  <p className="text-xs text-slate-500 mb-1">Blade</p>
                  <p className="text-sm font-semibold text-slate-800">
                    {getBladeName(selectedModel.blade_id)}
                  </p>
                </div>
                <div className="rounded-lg border border-slate-200 p-4">
                  <p className="text-xs text-slate-500 mb-1">Controller</p>
                  <p className="text-sm font-semibold text-slate-800">
                    {getControllerName(selectedModel.controller_id)}
                  </p>
                </div>
              </div>

              {/* Parameters */}
              <div className="grid grid-cols-2 gap-4">
                {[
                  ['Number of Blades', selectedModel.num_blades],
                  ['Hub Mass', `${selectedModel.hub_mass_kg.toLocaleString()} kg`],
                  ['Nacelle Mass', `${selectedModel.nacelle_mass_kg.toLocaleString()} kg`],
                  ['Shaft Tilt', `${selectedModel.shaft_tilt_deg} deg`],
                  ['Precone', `${selectedModel.precone_deg} deg`],
                  ['Overhang', `${selectedModel.overhang_m} m`],
                  ['Tower-to-Shaft', `${selectedModel.twr2shft_m} m`],
                  ['Gearbox Ratio', selectedModel.gearbox_ratio],
                  ['Drivetrain Eff.', `${(selectedModel.drivetrain_efficiency * 100).toFixed(1)}%`],
                  ['Gen Inertia', `${selectedModel.gen_inertia_kg_m2.toLocaleString()} kg*m2`],
                ].map(([label, value]) => (
                  <div
                    key={String(label)}
                    className="flex items-center justify-between rounded-md bg-slate-50 px-3 py-2"
                  >
                    <span className="text-xs text-slate-500">{label}</span>
                    <span className="text-sm font-medium text-slate-800 font-mono">
                      {value}
                    </span>
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
