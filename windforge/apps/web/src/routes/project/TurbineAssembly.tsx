import { useEffect, useState, useCallback } from 'react';
import { useParams } from 'react-router-dom';
import { turbineModelsApi, towersApi, bladesApi, controllersApi } from '@/api/client';
import type { TurbineModel, TurbineModelCreate, Tower, Blade, Controller } from '@/types';
import { Boxes, Plus, Loader2, Trash2, Save, X } from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';

const DEFAULT_DOF_FLAGS: Record<string, boolean> = {
  FlapDOF1: true,
  FlapDOF2: true,
  EdgeDOF: true,
  DrTrDOF: true,
  GenDOF: true,
  TwFADOF1: true,
  TwFADOF2: true,
  TwSSDOF1: true,
  TwSSDOF2: true,
  YawDOF: false,
};

interface EditableModel {
  name: string;
  tower_id: string;
  blade_id: string;
  controller_id: string;
  gearbox_ratio: number;
  generator_inertia: number;
  drivetrain_stiffness: number;
  drivetrain_damping: number;
  hub_mass: number;
  hub_inertia: number;
  nacelle_mass: number;
  nacelle_inertia: number;
  overhang: number;
  shaft_tilt: number;
  precone: number;
  rotor_speed_rated: number;
  dof_flags: Record<string, boolean>;
}

function modelToEditable(model: TurbineModel): EditableModel {
  return {
    name: model.name,
    tower_id: model.tower_id || '',
    blade_id: model.blade_id || '',
    controller_id: model.controller_id || '',
    gearbox_ratio: model.gearbox_ratio ?? 0,
    generator_inertia: model.generator_inertia ?? 0,
    drivetrain_stiffness: model.drivetrain_stiffness ?? 0,
    drivetrain_damping: model.drivetrain_damping ?? 0,
    hub_mass: model.hub_mass ?? 0,
    hub_inertia: model.hub_inertia ?? 0,
    nacelle_mass: model.nacelle_mass ?? 0,
    nacelle_inertia: model.nacelle_inertia ?? 0,
    overhang: model.overhang ?? 0,
    shaft_tilt: model.shaft_tilt ?? 0,
    precone: model.precone ?? 0,
    rotor_speed_rated: model.rotor_speed_rated ?? 0,
    dof_flags: model.dof_flags ? { ...model.dof_flags } : { ...DEFAULT_DOF_FLAGS },
  };
}

export default function TurbineAssembly() {
  const { projectId } = useParams<{ projectId: string }>();
  const [models, setModels] = useState<TurbineModel[]>([]);
  const [towers, setTowers] = useState<Tower[]>([]);
  const [blades, setBlades] = useState<Blade[]>([]);
  const [controllers, setControllers] = useState<Controller[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedModel, setSelectedModel] = useState<TurbineModel | null>(null);
  const [editForm, setEditForm] = useState<EditableModel | null>(null);
  const [isSaving, setIsSaving] = useState(false);
  const [showCreateModal, setShowCreateModal] = useState(false);
  const [newName, setNewName] = useState('');
  const [isCreating, setIsCreating] = useState(false);

  const loadData = useCallback(async () => {
    if (!projectId) return;
    try {
      const [modelsData, towersData, bladesData, controllersData] = await Promise.all([
        turbineModelsApi.list(projectId),
        towersApi.list(projectId),
        bladesApi.list(projectId),
        controllersApi.list(projectId),
      ]);
      setModels(modelsData);
      setTowers(towersData);
      setBlades(bladesData);
      setControllers(controllersData);
      return modelsData;
    } catch (err) {
      toast.error('Failed to load data');
      return [];
    }
  }, [projectId]);

  useEffect(() => {
    setIsLoading(true);
    loadData()
      .then((data) => {
        if (data && data.length > 0) {
          setSelectedModel(data[0]);
          setEditForm(modelToEditable(data[0]));
        }
      })
      .finally(() => setIsLoading(false));
  }, [loadData]);

  const handleSelect = (model: TurbineModel) => {
    setSelectedModel(model);
    setEditForm(modelToEditable(model));
  };

  const handleCreate = async () => {
    if (!projectId || !newName.trim()) return;
    setIsCreating(true);
    try {
      const created = await turbineModelsApi.create(projectId, {
        name: newName.trim(),
      } as TurbineModelCreate);
      toast.success('Turbine model created');
      setShowCreateModal(false);
      setNewName('');
      const data = await loadData();
      if (data) {
        const found = data.find((m) => m.id === created.id);
        if (found) {
          setSelectedModel(found);
          setEditForm(modelToEditable(found));
        }
      }
    } catch (err) {
      toast.error('Failed to create turbine model');
    } finally {
      setIsCreating(false);
    }
  };

  const handleSave = async () => {
    if (!projectId || !selectedModel || !editForm) return;
    setIsSaving(true);
    try {
      const payload: Partial<TurbineModelCreate> = {
        name: editForm.name,
        tower_id: editForm.tower_id || undefined,
        blade_id: editForm.blade_id || undefined,
        controller_id: editForm.controller_id || undefined,
        gearbox_ratio: editForm.gearbox_ratio,
        generator_inertia: editForm.generator_inertia,
        drivetrain_stiffness: editForm.drivetrain_stiffness,
        drivetrain_damping: editForm.drivetrain_damping,
        hub_mass: editForm.hub_mass,
        hub_inertia: editForm.hub_inertia,
        nacelle_mass: editForm.nacelle_mass,
        nacelle_inertia: editForm.nacelle_inertia,
        overhang: editForm.overhang,
        shaft_tilt: editForm.shaft_tilt,
        precone: editForm.precone,
        rotor_speed_rated: editForm.rotor_speed_rated,
        dof_flags: editForm.dof_flags,
      };
      const updated = await turbineModelsApi.update(projectId, selectedModel.id, payload);
      toast.success('Turbine model saved');
      setSelectedModel(updated);
      setEditForm(modelToEditable(updated));
      await loadData();
    } catch (err) {
      toast.error('Failed to save turbine model');
    } finally {
      setIsSaving(false);
    }
  };

  const handleDelete = async () => {
    if (!projectId || !selectedModel) return;
    if (!window.confirm(`Delete turbine model "${selectedModel.name}"?`)) return;
    try {
      await turbineModelsApi.delete(projectId, selectedModel.id);
      toast.success('Turbine model deleted');
      setSelectedModel(null);
      setEditForm(null);
      await loadData();
    } catch (err) {
      toast.error('Failed to delete turbine model');
    }
  };

  const updateField = (field: keyof EditableModel, value: any) => {
    if (!editForm) return;
    setEditForm({ ...editForm, [field]: value });
  };

  const toggleDof = (key: string) => {
    if (!editForm) return;
    setEditForm({
      ...editForm,
      dof_flags: { ...editForm.dof_flags, [key]: !editForm.dof_flags[key] },
    });
  };

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
          <h2 className="text-xl font-bold text-slate-800">Turbine Assembly</h2>
          <p className="text-sm text-slate-500">
            Combine tower, blade, and controller into a complete turbine model
          </p>
        </div>
        <button
          onClick={() => setShowCreateModal(true)}
          className="btn-primary flex items-center gap-2"
        >
          <Plus className="h-4 w-4" />
          New Turbine Model
        </button>
      </div>

      {/* Create Modal */}
      {showCreateModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40">
          <div className="bg-white rounded-xl shadow-xl p-6 w-full max-w-md">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-lg font-semibold text-slate-800">New Turbine Model</h3>
              <button onClick={() => setShowCreateModal(false)} className="text-slate-400 hover:text-slate-600">
                <X className="h-5 w-5" />
              </button>
            </div>
            <div>
              <label className="block text-sm font-medium text-slate-700 mb-1">Name</label>
              <input
                type="text"
                value={newName}
                onChange={(e) => setNewName(e.target.value)}
                placeholder="e.g., NREL 5MW"
                className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
              />
            </div>
            <div className="flex justify-end gap-3 mt-6">
              <button
                onClick={() => setShowCreateModal(false)}
                className="px-4 py-2 text-sm font-medium text-slate-700 bg-slate-100 rounded-lg hover:bg-slate-200"
              >
                Cancel
              </button>
              <button
                onClick={handleCreate}
                disabled={isCreating || !newName.trim()}
                className="px-4 py-2 text-sm font-medium text-white bg-accent-600 rounded-lg hover:bg-accent-700 disabled:opacity-50 flex items-center gap-2"
              >
                {isCreating && <Loader2 className="h-4 w-4 animate-spin" />}
                Create
              </button>
            </div>
          </div>
        </div>
      )}

      {models.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Boxes className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">No turbine assemblies</h3>
          <p className="mt-1 text-sm text-slate-500 max-w-sm text-center">
            Define tower, blade, and controller first, then assemble them into a complete turbine model.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Model list */}
          <div className="space-y-2">
            {models.map((model) => (
              <button
                key={model.id}
                onClick={() => handleSelect(model)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedModel?.id === model.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <span className="font-medium text-slate-800">{model.name}</span>
                <div className="mt-1 text-xs text-slate-500">
                  GR {model.gearbox_ratio ?? '--'}
                </div>
              </button>
            ))}
          </div>

          {/* Model detail */}
          {selectedModel && editForm && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6 space-y-6">
              <div className="flex items-center justify-between">
                <h3 className="text-lg font-semibold text-slate-800">Edit Turbine Model</h3>
                <div className="flex items-center gap-2">
                  <button
                    onClick={handleSave}
                    disabled={isSaving}
                    className="flex items-center gap-2 px-4 py-2 text-sm font-medium text-white bg-accent-600 rounded-lg hover:bg-accent-700 disabled:opacity-50"
                  >
                    {isSaving ? <Loader2 className="h-4 w-4 animate-spin" /> : <Save className="h-4 w-4" />}
                    Save
                  </button>
                  <button
                    onClick={handleDelete}
                    className="flex items-center gap-2 px-4 py-2 text-sm font-medium text-danger-600 bg-danger-50 rounded-lg hover:bg-danger-100"
                  >
                    <Trash2 className="h-4 w-4" />
                    Delete
                  </button>
                </div>
              </div>

              {/* Name */}
              <div>
                <label className="block text-xs font-medium text-slate-500 mb-1">Name</label>
                <input
                  type="text"
                  value={editForm.name}
                  onChange={(e) => updateField('name', e.target.value)}
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                />
              </div>

              {/* Component dropdowns */}
              <div>
                <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider mb-3">Components</h4>
                <div className="grid grid-cols-3 gap-4">
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Tower</label>
                    <select
                      value={editForm.tower_id}
                      onChange={(e) => updateField('tower_id', e.target.value)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                    >
                      <option value="">-- Select --</option>
                      {towers.map((t) => (
                        <option key={t.id} value={t.id}>{t.name}</option>
                      ))}
                    </select>
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Blade</label>
                    <select
                      value={editForm.blade_id}
                      onChange={(e) => updateField('blade_id', e.target.value)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                    >
                      <option value="">-- Select --</option>
                      {blades.map((b) => (
                        <option key={b.id} value={b.id}>{b.name}</option>
                      ))}
                    </select>
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Controller</label>
                    <select
                      value={editForm.controller_id}
                      onChange={(e) => updateField('controller_id', e.target.value)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                    >
                      <option value="">-- Select --</option>
                      {controllers.map((c) => (
                        <option key={c.id} value={c.id}>{c.name}</option>
                      ))}
                    </select>
                  </div>
                </div>
              </div>

              {/* Drivetrain */}
              <div>
                <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider mb-3">Drivetrain</h4>
                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Gearbox Ratio</label>
                    <input type="number" step="0.01" value={editForm.gearbox_ratio}
                      onChange={(e) => updateField('gearbox_ratio', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Generator Inertia (kg*m2)</label>
                    <input type="number" step="0.1" value={editForm.generator_inertia}
                      onChange={(e) => updateField('generator_inertia', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Drivetrain Stiffness (Nm/rad)</label>
                    <input type="number" step="1e6" value={editForm.drivetrain_stiffness}
                      onChange={(e) => updateField('drivetrain_stiffness', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Drivetrain Damping (Nm/(rad/s))</label>
                    <input type="number" step="1e4" value={editForm.drivetrain_damping}
                      onChange={(e) => updateField('drivetrain_damping', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                </div>
              </div>

              {/* Hub / Nacelle */}
              <div>
                <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider mb-3">Hub / Nacelle</h4>
                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Hub Mass (kg)</label>
                    <input type="number" step="1" value={editForm.hub_mass}
                      onChange={(e) => updateField('hub_mass', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Hub Inertia (kg*m2)</label>
                    <input type="number" step="1" value={editForm.hub_inertia}
                      onChange={(e) => updateField('hub_inertia', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Nacelle Mass (kg)</label>
                    <input type="number" step="1" value={editForm.nacelle_mass}
                      onChange={(e) => updateField('nacelle_mass', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Nacelle Inertia (kg*m2)</label>
                    <input type="number" step="1" value={editForm.nacelle_inertia}
                      onChange={(e) => updateField('nacelle_inertia', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Overhang (m)</label>
                    <input type="number" step="0.01" value={editForm.overhang}
                      onChange={(e) => updateField('overhang', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Shaft Tilt (deg)</label>
                    <input type="number" step="0.1" value={editForm.shaft_tilt}
                      onChange={(e) => updateField('shaft_tilt', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Precone (deg)</label>
                    <input type="number" step="0.1" value={editForm.precone}
                      onChange={(e) => updateField('precone', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs font-medium text-slate-500 mb-1">Rotor Speed Rated (rpm)</label>
                    <input type="number" step="0.1" value={editForm.rotor_speed_rated}
                      onChange={(e) => updateField('rotor_speed_rated', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                </div>
              </div>

              {/* DOF Flags */}
              <div>
                <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider mb-3">
                  Degrees of Freedom
                </h4>
                <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
                  {Object.entries(editForm.dof_flags).map(([key, val]) => (
                    <label
                      key={key}
                      className="flex items-center gap-2 rounded-lg bg-slate-50 px-3 py-2 cursor-pointer hover:bg-slate-100"
                    >
                      <input
                        type="checkbox"
                        checked={val}
                        onChange={() => toggleDof(key)}
                        className="rounded border-slate-300 text-accent-600 focus:ring-accent-500"
                      />
                      <span className="text-sm text-slate-700">{key}</span>
                    </label>
                  ))}
                </div>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
