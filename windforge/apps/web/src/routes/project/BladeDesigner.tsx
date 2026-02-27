import { useEffect, useState, useCallback } from 'react';
import { useParams } from 'react-router-dom';
import { bladesApi } from '@/api/client';
import type { Blade, BladeCreate, BladeStructuralStation, BladeAeroStation } from '@/types';
import { Fan, Plus, Loader2, Trash2, Save, X } from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';

interface EditableBlade {
  name: string;
  blade_length: number;
  blade_flap_damping: number;
  blade_edge_damping: number;
  structural_stations: BladeStructuralStation[];
  aero_stations: BladeAeroStation[];
}

function bladeToEditable(blade: Blade): EditableBlade {
  return {
    name: blade.name,
    blade_length: blade.blade_length,
    blade_flap_damping: blade.blade_flap_damping,
    blade_edge_damping: blade.blade_edge_damping,
    structural_stations: blade.structural_stations
      ? blade.structural_stations.map((s) => ({ ...s }))
      : [],
    aero_stations: blade.aero_stations
      ? blade.aero_stations.map((s) => ({ ...s }))
      : [],
  };
}

export default function BladeDesigner() {
  const { projectId } = useParams<{ projectId: string }>();
  const [blades, setBlades] = useState<Blade[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedBlade, setSelectedBlade] = useState<Blade | null>(null);
  const [editForm, setEditForm] = useState<EditableBlade | null>(null);
  const [activeTab, setActiveTab] = useState<'structural' | 'aero'>('structural');
  const [isSaving, setIsSaving] = useState(false);
  const [showCreateModal, setShowCreateModal] = useState(false);
  const [newName, setNewName] = useState('');
  const [newLength, setNewLength] = useState<number>(61.5);
  const [isCreating, setIsCreating] = useState(false);

  const loadBlades = useCallback(async () => {
    if (!projectId) return;
    try {
      const data = await bladesApi.list(projectId);
      setBlades(data);
      return data;
    } catch (err) {
      toast.error('Failed to load blades');
      return [];
    }
  }, [projectId]);

  useEffect(() => {
    setIsLoading(true);
    loadBlades()
      .then((data) => {
        if (data && data.length > 0) {
          setSelectedBlade(data[0]);
          setEditForm(bladeToEditable(data[0]));
        }
      })
      .finally(() => setIsLoading(false));
  }, [loadBlades]);

  const handleSelect = (blade: Blade) => {
    setSelectedBlade(blade);
    setEditForm(bladeToEditable(blade));
  };

  const handleCreate = async () => {
    if (!projectId || !newName.trim()) return;
    setIsCreating(true);
    try {
      const created = await bladesApi.create(projectId, {
        name: newName.trim(),
        blade_length: newLength,
      } as BladeCreate);
      toast.success('Blade created');
      setShowCreateModal(false);
      setNewName('');
      setNewLength(61.5);
      const data = await loadBlades();
      if (data) {
        const found = data.find((b) => b.id === created.id);
        if (found) {
          setSelectedBlade(found);
          setEditForm(bladeToEditable(found));
        }
      }
    } catch (err) {
      toast.error('Failed to create blade');
    } finally {
      setIsCreating(false);
    }
  };

  const handleSave = async () => {
    if (!projectId || !selectedBlade || !editForm) return;
    setIsSaving(true);
    try {
      const payload: Partial<BladeCreate> = {
        name: editForm.name,
        blade_length: editForm.blade_length,
        blade_flap_damping: editForm.blade_flap_damping,
        blade_edge_damping: editForm.blade_edge_damping,
        structural_stations: editForm.structural_stations,
        aero_stations: editForm.aero_stations,
      };
      const updated = await bladesApi.update(projectId, selectedBlade.id, payload);
      toast.success('Blade saved');
      setSelectedBlade(updated);
      setEditForm(bladeToEditable(updated));
      await loadBlades();
    } catch (err) {
      toast.error('Failed to save blade');
    } finally {
      setIsSaving(false);
    }
  };

  const handleDelete = async () => {
    if (!projectId || !selectedBlade) return;
    if (!window.confirm(`Delete blade "${selectedBlade.name}"?`)) return;
    try {
      await bladesApi.delete(projectId, selectedBlade.id);
      toast.success('Blade deleted');
      setSelectedBlade(null);
      setEditForm(null);
      await loadBlades();
    } catch (err) {
      toast.error('Failed to delete blade');
    }
  };

  const updateField = (field: keyof EditableBlade, value: any) => {
    if (!editForm) return;
    setEditForm({ ...editForm, [field]: value });
  };

  const updateStructStation = (index: number, field: keyof BladeStructuralStation, value: number) => {
    if (!editForm) return;
    const stations = [...editForm.structural_stations];
    stations[index] = { ...stations[index], [field]: value };
    setEditForm({ ...editForm, structural_stations: stations });
  };

  const updateAeroStation = (index: number, field: keyof BladeAeroStation, value: any) => {
    if (!editForm) return;
    const stations = [...editForm.aero_stations];
    stations[index] = { ...stations[index], [field]: value };
    setEditForm({ ...editForm, aero_stations: stations });
  };

  const addStructStation = () => {
    if (!editForm) return;
    const newStation: BladeStructuralStation = {
      frac: editForm.structural_stations.length > 0 ? 1.0 : 0.0,
      pitch_axis: 0.5,
      struct_twist: 0,
      mass_den: 0,
      flap_stiff: 0,
      edge_stiff: 0,
    };
    setEditForm({ ...editForm, structural_stations: [...editForm.structural_stations, newStation] });
  };

  const removeStructStation = (index: number) => {
    if (!editForm) return;
    setEditForm({
      ...editForm,
      structural_stations: editForm.structural_stations.filter((_, i) => i !== index),
    });
  };

  const addAeroStation = () => {
    if (!editForm) return;
    const newStation: BladeAeroStation = {
      frac: editForm.aero_stations.length > 0 ? 1.0 : 0.0,
      chord: 0,
      aero_twist: 0,
      airfoil_id: '',
      aero_center: 0.25,
    };
    setEditForm({ ...editForm, aero_stations: [...editForm.aero_stations, newStation] });
  };

  const removeAeroStation = (index: number) => {
    if (!editForm) return;
    setEditForm({
      ...editForm,
      aero_stations: editForm.aero_stations.filter((_, i) => i !== index),
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
          <h2 className="text-xl font-bold text-slate-800">Blade Designer</h2>
          <p className="text-sm text-slate-500">
            Define blade structural and aerodynamic properties
          </p>
        </div>
        <button
          onClick={() => setShowCreateModal(true)}
          className="btn-primary flex items-center gap-2"
        >
          <Plus className="h-4 w-4" />
          New Blade
        </button>
      </div>

      {/* Create Modal */}
      {showCreateModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40">
          <div className="bg-white rounded-xl shadow-xl p-6 w-full max-w-md">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-lg font-semibold text-slate-800">New Blade</h3>
              <button onClick={() => setShowCreateModal(false)} className="text-slate-400 hover:text-slate-600">
                <X className="h-5 w-5" />
              </button>
            </div>
            <div className="space-y-4">
              <div>
                <label className="block text-sm font-medium text-slate-700 mb-1">Name</label>
                <input
                  type="text"
                  value={newName}
                  onChange={(e) => setNewName(e.target.value)}
                  placeholder="e.g., NREL 5MW Blade"
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                />
              </div>
              <div>
                <label className="block text-sm font-medium text-slate-700 mb-1">Blade Length (m)</label>
                <input
                  type="number"
                  value={newLength}
                  onChange={(e) => setNewLength(parseFloat(e.target.value) || 0)}
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                />
              </div>
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

      {blades.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Fan className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">No blades defined</h3>
          <p className="mt-1 text-sm text-slate-500">
            Create your first blade design to define structural and aerodynamic properties.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Blade list */}
          <div className="space-y-2">
            {blades.map((blade) => (
              <button
                key={blade.id}
                onClick={() => handleSelect(blade)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedBlade?.id === blade.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <div className="flex items-center justify-between">
                  <span className="font-medium text-slate-800">{blade.name}</span>
                </div>
                <div className="mt-1 text-xs text-slate-500">
                  {blade.blade_length}m long
                </div>
              </button>
            ))}
          </div>

          {/* Blade detail */}
          {selectedBlade && editForm && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-6">
                <h3 className="text-lg font-semibold text-slate-800">Edit Blade</h3>
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

              {/* Basic fields */}
              <div className="grid grid-cols-2 gap-4 mb-6">
                <div className="col-span-2">
                  <label className="block text-xs font-medium text-slate-500 mb-1">Name</label>
                  <input
                    type="text"
                    value={editForm.name}
                    onChange={(e) => updateField('name', e.target.value)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Blade Length (m)</label>
                  <input
                    type="number"
                    value={editForm.blade_length}
                    onChange={(e) => updateField('blade_length', parseFloat(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Flap Damping (%)</label>
                  <input
                    type="number"
                    step="0.01"
                    value={editForm.blade_flap_damping}
                    onChange={(e) => updateField('blade_flap_damping', parseFloat(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Edge Damping (%)</label>
                  <input
                    type="number"
                    step="0.01"
                    value={editForm.blade_edge_damping}
                    onChange={(e) => updateField('blade_edge_damping', parseFloat(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
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
                  Structural ({editForm.structural_stations.length})
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
                  Aerodynamic ({editForm.aero_stations.length})
                </button>
              </div>

              {/* Structural tab */}
              {activeTab === 'structural' && (
                <>
                  <div className="flex justify-end mb-2">
                    <button
                      onClick={addStructStation}
                      className="flex items-center gap-1 px-3 py-1.5 text-xs font-medium text-accent-700 bg-accent-50 rounded-lg hover:bg-accent-100"
                    >
                      <Plus className="h-3.5 w-3.5" />
                      Add Station
                    </button>
                  </div>
                  {editForm.structural_stations.length > 0 ? (
                    <div className="overflow-x-auto rounded-lg border border-slate-200">
                      <table className="min-w-full text-sm">
                        <thead className="bg-slate-50">
                          <tr>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Frac</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Pitch Axis</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Struct Twist (deg)</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Mass Den (kg/m)</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Flap Stiff (Nm2)</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Edge Stiff (Nm2)</th>
                            <th className="px-3 py-2 w-10"></th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-slate-100">
                          {editForm.structural_stations.map((s, i) => (
                            <tr key={i} className="hover:bg-slate-50">
                              <td className="px-2 py-1">
                                <input type="number" step="0.001" value={s.frac}
                                  onChange={(e) => updateStructStation(i, 'frac', parseFloat(e.target.value) || 0)}
                                  className="w-20 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <input type="number" step="0.01" value={s.pitch_axis}
                                  onChange={(e) => updateStructStation(i, 'pitch_axis', parseFloat(e.target.value) || 0)}
                                  className="w-20 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <input type="number" step="0.1" value={s.struct_twist}
                                  onChange={(e) => updateStructStation(i, 'struct_twist', parseFloat(e.target.value) || 0)}
                                  className="w-24 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <input type="number" step="0.1" value={s.mass_den}
                                  onChange={(e) => updateStructStation(i, 'mass_den', parseFloat(e.target.value) || 0)}
                                  className="w-24 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <input type="number" step="1e6" value={s.flap_stiff}
                                  onChange={(e) => updateStructStation(i, 'flap_stiff', parseFloat(e.target.value) || 0)}
                                  className="w-28 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <input type="number" step="1e6" value={s.edge_stiff}
                                  onChange={(e) => updateStructStation(i, 'edge_stiff', parseFloat(e.target.value) || 0)}
                                  className="w-28 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <button onClick={() => removeStructStation(i)} className="text-slate-400 hover:text-danger-500">
                                  <Trash2 className="h-3.5 w-3.5" />
                                </button>
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  ) : (
                    <div className="rounded-lg border-2 border-dashed border-slate-200 p-8 text-center">
                      <p className="text-sm text-slate-500">No structural stations. Click "Add Station" to begin.</p>
                    </div>
                  )}
                </>
              )}

              {/* Aero tab */}
              {activeTab === 'aero' && (
                <>
                  <div className="flex justify-end mb-2">
                    <button
                      onClick={addAeroStation}
                      className="flex items-center gap-1 px-3 py-1.5 text-xs font-medium text-accent-700 bg-accent-50 rounded-lg hover:bg-accent-100"
                    >
                      <Plus className="h-3.5 w-3.5" />
                      Add Station
                    </button>
                  </div>
                  {editForm.aero_stations.length > 0 ? (
                    <div className="overflow-x-auto rounded-lg border border-slate-200">
                      <table className="min-w-full text-sm">
                        <thead className="bg-slate-50">
                          <tr>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Frac</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Chord (m)</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Aero Twist (deg)</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Airfoil ID</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Aero Center</th>
                            <th className="px-3 py-2 w-10"></th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-slate-100">
                          {editForm.aero_stations.map((s, i) => (
                            <tr key={i} className="hover:bg-slate-50">
                              <td className="px-2 py-1">
                                <input type="number" step="0.001" value={s.frac}
                                  onChange={(e) => updateAeroStation(i, 'frac', parseFloat(e.target.value) || 0)}
                                  className="w-20 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <input type="number" step="0.001" value={s.chord}
                                  onChange={(e) => updateAeroStation(i, 'chord', parseFloat(e.target.value) || 0)}
                                  className="w-24 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <input type="number" step="0.01" value={s.aero_twist}
                                  onChange={(e) => updateAeroStation(i, 'aero_twist', parseFloat(e.target.value) || 0)}
                                  className="w-24 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <input type="text" value={s.airfoil_id}
                                  onChange={(e) => updateAeroStation(i, 'airfoil_id', e.target.value)}
                                  className="w-32 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <input type="number" step="0.01" value={s.aero_center}
                                  onChange={(e) => updateAeroStation(i, 'aero_center', parseFloat(e.target.value) || 0)}
                                  className="w-20 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                              </td>
                              <td className="px-2 py-1">
                                <button onClick={() => removeAeroStation(i)} className="text-slate-400 hover:text-danger-500">
                                  <Trash2 className="h-3.5 w-3.5" />
                                </button>
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  ) : (
                    <div className="rounded-lg border-2 border-dashed border-slate-200 p-8 text-center">
                      <p className="text-sm text-slate-500">No aero stations. Click "Add Station" to begin.</p>
                    </div>
                  )}
                </>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
