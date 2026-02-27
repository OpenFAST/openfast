import { useEffect, useState, useCallback } from 'react';
import { useParams } from 'react-router-dom';
import { towersApi } from '@/api/client';
import type { Tower, TowerCreate, TowerStation } from '@/types';
import { Building2, Plus, Loader2, Trash2, Save, X } from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';

interface EditableTower {
  name: string;
  tower_height: number;
  tower_base_height: number;
  tower_fa_damping_1: number;
  tower_fa_damping_2: number;
  tower_ss_damping_1: number;
  tower_ss_damping_2: number;
  stations: TowerStation[];
}

function towerToEditable(tower: Tower): EditableTower {
  return {
    name: tower.name,
    tower_height: tower.tower_height,
    tower_base_height: tower.tower_base_height,
    tower_fa_damping_1: tower.tower_fa_damping_1,
    tower_fa_damping_2: tower.tower_fa_damping_2,
    tower_ss_damping_1: tower.tower_ss_damping_1,
    tower_ss_damping_2: tower.tower_ss_damping_2,
    stations: tower.stations ? [...tower.stations.map((s) => ({ ...s }))] : [],
  };
}

export default function TowerDesigner() {
  const { projectId } = useParams<{ projectId: string }>();
  const [towers, setTowers] = useState<Tower[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedTower, setSelectedTower] = useState<Tower | null>(null);
  const [editForm, setEditForm] = useState<EditableTower | null>(null);
  const [isSaving, setIsSaving] = useState(false);
  const [showCreateModal, setShowCreateModal] = useState(false);
  const [newName, setNewName] = useState('');
  const [newHeight, setNewHeight] = useState<number>(80);
  const [isCreating, setIsCreating] = useState(false);

  const loadTowers = useCallback(async () => {
    if (!projectId) return;
    try {
      const data = await towersApi.list(projectId);
      setTowers(data);
      return data;
    } catch (err) {
      toast.error('Failed to load towers');
      return [];
    }
  }, [projectId]);

  useEffect(() => {
    setIsLoading(true);
    loadTowers()
      .then((data) => {
        if (data && data.length > 0) {
          setSelectedTower(data[0]);
          setEditForm(towerToEditable(data[0]));
        }
      })
      .finally(() => setIsLoading(false));
  }, [loadTowers]);

  const handleSelect = (tower: Tower) => {
    setSelectedTower(tower);
    setEditForm(towerToEditable(tower));
  };

  const handleCreate = async () => {
    if (!projectId || !newName.trim()) return;
    setIsCreating(true);
    try {
      const created = await towersApi.create(projectId, {
        name: newName.trim(),
        tower_height: newHeight,
      });
      toast.success('Tower created');
      setShowCreateModal(false);
      setNewName('');
      setNewHeight(80);
      const data = await loadTowers();
      if (data) {
        const found = data.find((t) => t.id === created.id);
        if (found) {
          setSelectedTower(found);
          setEditForm(towerToEditable(found));
        }
      }
    } catch (err) {
      toast.error('Failed to create tower');
    } finally {
      setIsCreating(false);
    }
  };

  const handleSave = async () => {
    if (!projectId || !selectedTower || !editForm) return;
    setIsSaving(true);
    try {
      const payload: Partial<TowerCreate> = {
        name: editForm.name,
        tower_height: editForm.tower_height,
        tower_base_height: editForm.tower_base_height,
        tower_fa_damping_1: editForm.tower_fa_damping_1,
        tower_fa_damping_2: editForm.tower_fa_damping_2,
        tower_ss_damping_1: editForm.tower_ss_damping_1,
        tower_ss_damping_2: editForm.tower_ss_damping_2,
        stations: editForm.stations,
      };
      const updated = await towersApi.update(projectId, selectedTower.id, payload);
      toast.success('Tower saved');
      setSelectedTower(updated);
      setEditForm(towerToEditable(updated));
      await loadTowers();
    } catch (err) {
      toast.error('Failed to save tower');
    } finally {
      setIsSaving(false);
    }
  };

  const handleDelete = async () => {
    if (!projectId || !selectedTower) return;
    if (!window.confirm(`Delete tower "${selectedTower.name}"?`)) return;
    try {
      await towersApi.delete(projectId, selectedTower.id);
      toast.success('Tower deleted');
      setSelectedTower(null);
      setEditForm(null);
      await loadTowers();
    } catch (err) {
      toast.error('Failed to delete tower');
    }
  };

  const updateField = (field: keyof EditableTower, value: any) => {
    if (!editForm) return;
    setEditForm({ ...editForm, [field]: value });
  };

  const updateStation = (index: number, field: keyof TowerStation, value: number) => {
    if (!editForm) return;
    const stations = [...editForm.stations];
    stations[index] = { ...stations[index], [field]: value };
    setEditForm({ ...editForm, stations });
  };

  const addStation = () => {
    if (!editForm) return;
    const newStation: TowerStation = {
      frac: editForm.stations.length > 0 ? 1.0 : 0.0,
      mass_den: 0,
      fa_stiff: 0,
      ss_stiff: 0,
      outer_diameter: 0,
      wall_thickness: 0,
    };
    setEditForm({ ...editForm, stations: [...editForm.stations, newStation] });
  };

  const removeStation = (index: number) => {
    if (!editForm) return;
    const stations = editForm.stations.filter((_, i) => i !== index);
    setEditForm({ ...editForm, stations });
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
          <h2 className="text-xl font-bold text-slate-800">Tower Designer</h2>
          <p className="text-sm text-slate-500">
            Define tower geometry, materials, and structural stations
          </p>
        </div>
        <button
          onClick={() => setShowCreateModal(true)}
          className="btn-primary flex items-center gap-2"
        >
          <Plus className="h-4 w-4" />
          New Tower
        </button>
      </div>

      {/* Create Modal */}
      {showCreateModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40">
          <div className="bg-white rounded-xl shadow-xl p-6 w-full max-w-md">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-lg font-semibold text-slate-800">New Tower</h3>
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
                  placeholder="e.g., NREL 5MW Tower"
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                />
              </div>
              <div>
                <label className="block text-sm font-medium text-slate-700 mb-1">Tower Height (m)</label>
                <input
                  type="number"
                  value={newHeight}
                  onChange={(e) => setNewHeight(parseFloat(e.target.value) || 0)}
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

      {towers.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Building2 className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">No towers defined</h3>
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
                onClick={() => handleSelect(tower)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedTower?.id === tower.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <div className="flex items-center justify-between">
                  <span className="font-medium text-slate-800">{tower.name}</span>
                </div>
                <div className="mt-1 text-xs text-slate-500">
                  {tower.tower_height}m tall &middot;{' '}
                  {tower.stations ? tower.stations.length : 0} stations
                </div>
              </button>
            ))}
          </div>

          {/* Tower detail */}
          {selectedTower && editForm && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-6">
                <h3 className="text-lg font-semibold text-slate-800">Edit Tower</h3>
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
                  <label className="block text-xs font-medium text-slate-500 mb-1">Tower Height (m)</label>
                  <input
                    type="number"
                    value={editForm.tower_height}
                    onChange={(e) => updateField('tower_height', parseFloat(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Base Height (m)</label>
                  <input
                    type="number"
                    value={editForm.tower_base_height}
                    onChange={(e) => updateField('tower_base_height', parseFloat(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
              </div>

              {/* Damping */}
              <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider mb-3">
                Damping Coefficients
              </h4>
              <div className="grid grid-cols-2 gap-4 mb-6">
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">FA Damping 1 (%)</label>
                  <input
                    type="number"
                    step="0.01"
                    value={editForm.tower_fa_damping_1}
                    onChange={(e) => updateField('tower_fa_damping_1', parseFloat(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">FA Damping 2 (%)</label>
                  <input
                    type="number"
                    step="0.01"
                    value={editForm.tower_fa_damping_2}
                    onChange={(e) => updateField('tower_fa_damping_2', parseFloat(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">SS Damping 1 (%)</label>
                  <input
                    type="number"
                    step="0.01"
                    value={editForm.tower_ss_damping_1}
                    onChange={(e) => updateField('tower_ss_damping_1', parseFloat(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">SS Damping 2 (%)</label>
                  <input
                    type="number"
                    step="0.01"
                    value={editForm.tower_ss_damping_2}
                    onChange={(e) => updateField('tower_ss_damping_2', parseFloat(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
              </div>

              {/* Stations table */}
              <div className="flex items-center justify-between mb-3">
                <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider">
                  Structural Stations ({editForm.stations.length})
                </h4>
                <button
                  onClick={addStation}
                  className="flex items-center gap-1 px-3 py-1.5 text-xs font-medium text-accent-700 bg-accent-50 rounded-lg hover:bg-accent-100"
                >
                  <Plus className="h-3.5 w-3.5" />
                  Add Station
                </button>
              </div>

              {editForm.stations.length > 0 ? (
                <div className="overflow-x-auto rounded-lg border border-slate-200">
                  <table className="min-w-full text-sm">
                    <thead className="bg-slate-50">
                      <tr>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Frac</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Mass Den (kg/m)</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">FA Stiff (Nm2)</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">SS Stiff (Nm2)</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">OD (m)</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Wall Thick (m)</th>
                        <th className="px-3 py-2 w-10"></th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-100">
                      {editForm.stations.map((s, i) => (
                        <tr key={i} className="hover:bg-slate-50">
                          <td className="px-2 py-1">
                            <input
                              type="number"
                              step="0.001"
                              value={s.frac}
                              onChange={(e) => updateStation(i, 'frac', parseFloat(e.target.value) || 0)}
                              className="w-20 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </td>
                          <td className="px-2 py-1">
                            <input
                              type="number"
                              step="0.1"
                              value={s.mass_den}
                              onChange={(e) => updateStation(i, 'mass_den', parseFloat(e.target.value) || 0)}
                              className="w-24 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </td>
                          <td className="px-2 py-1">
                            <input
                              type="number"
                              step="1e6"
                              value={s.fa_stiff}
                              onChange={(e) => updateStation(i, 'fa_stiff', parseFloat(e.target.value) || 0)}
                              className="w-28 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </td>
                          <td className="px-2 py-1">
                            <input
                              type="number"
                              step="1e6"
                              value={s.ss_stiff}
                              onChange={(e) => updateStation(i, 'ss_stiff', parseFloat(e.target.value) || 0)}
                              className="w-28 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </td>
                          <td className="px-2 py-1">
                            <input
                              type="number"
                              step="0.01"
                              value={s.outer_diameter}
                              onChange={(e) => updateStation(i, 'outer_diameter', parseFloat(e.target.value) || 0)}
                              className="w-20 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </td>
                          <td className="px-2 py-1">
                            <input
                              type="number"
                              step="0.001"
                              value={s.wall_thickness}
                              onChange={(e) => updateStation(i, 'wall_thickness', parseFloat(e.target.value) || 0)}
                              className="w-20 rounded border border-slate-200 px-2 py-1 text-xs focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </td>
                          <td className="px-2 py-1">
                            <button
                              onClick={() => removeStation(i)}
                              className="text-slate-400 hover:text-danger-500"
                            >
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
                  <p className="text-sm text-slate-500">No stations defined. Click "Add Station" to begin.</p>
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
