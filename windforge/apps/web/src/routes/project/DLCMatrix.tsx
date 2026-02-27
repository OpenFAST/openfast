import { useEffect, useState, useCallback } from 'react';
import { useParams } from 'react-router-dom';
import { dlcDefinitionsApi, turbineModelsApi } from '@/api/client';
import type { DLCDefinition, DLCDefinitionCreate, DLCCaseSpec, TurbSimParams, TurbineModel } from '@/types';
import { Table2, Plus, Loader2, Trash2, Save, X } from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';

const DEFAULT_TURBSIM: TurbSimParams = {
  turbulence_model: 'IECKAI',
  iec_standard: '1-ED3',
  iec_turbc: 'B',
  grid_height: 150,
  grid_width: 150,
  num_grid_z: 25,
  num_grid_y: 25,
  time_step: 0.05,
  analysis_time: 630,
  ref_height: 90,
};

const DEFAULT_CASE: DLCCaseSpec = {
  dlc_number: '1.1',
  wind_speeds: [8, 10, 12, 14, 16],
  seeds: 6,
  yaw_misalignments: [0],
};

interface EditableDLC {
  name: string;
  turbine_model_id: string;
  dlc_cases: DLCCaseSpec[];
  turbsim_params: TurbSimParams;
}

function dlcToEditable(dlc: DLCDefinition): EditableDLC {
  return {
    name: dlc.name,
    turbine_model_id: dlc.turbine_model_id,
    dlc_cases: dlc.dlc_cases ? dlc.dlc_cases.map((c) => ({ ...c, wind_speeds: [...c.wind_speeds], yaw_misalignments: [...c.yaw_misalignments] })) : [],
    turbsim_params: dlc.turbsim_params ? { ...dlc.turbsim_params } : { ...DEFAULT_TURBSIM },
  };
}

function computeTotalCases(cases: DLCCaseSpec[]): number {
  return cases.reduce((sum, c) => sum + c.wind_speeds.length * c.seeds * c.yaw_misalignments.length, 0);
}

export default function DLCMatrix() {
  const { projectId } = useParams<{ projectId: string }>();
  const [dlcDefs, setDlcDefs] = useState<DLCDefinition[]>([]);
  const [turbineModels, setTurbineModels] = useState<TurbineModel[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedDlc, setSelectedDlc] = useState<DLCDefinition | null>(null);
  const [editForm, setEditForm] = useState<EditableDLC | null>(null);
  const [isSaving, setIsSaving] = useState(false);
  const [showCreateModal, setShowCreateModal] = useState(false);
  const [newName, setNewName] = useState('');
  const [newTurbineModelId, setNewTurbineModelId] = useState('');
  const [isCreating, setIsCreating] = useState(false);

  const loadData = useCallback(async () => {
    if (!projectId) return;
    try {
      const [defs, models] = await Promise.all([
        dlcDefinitionsApi.list(projectId),
        turbineModelsApi.list(projectId),
      ]);
      setDlcDefs(defs);
      setTurbineModels(models);
      return defs;
    } catch (err) {
      toast.error('Failed to load DLC definitions');
      return [];
    }
  }, [projectId]);

  useEffect(() => {
    setIsLoading(true);
    loadData()
      .then((data) => {
        if (data && data.length > 0) {
          setSelectedDlc(data[0]);
          setEditForm(dlcToEditable(data[0]));
        }
      })
      .finally(() => setIsLoading(false));
  }, [loadData]);

  const handleSelect = (dlc: DLCDefinition) => {
    setSelectedDlc(dlc);
    setEditForm(dlcToEditable(dlc));
  };

  const handleCreate = async () => {
    if (!projectId || !newName.trim() || !newTurbineModelId) return;
    setIsCreating(true);
    try {
      const created = await dlcDefinitionsApi.create(projectId, {
        name: newName.trim(),
        turbine_model_id: newTurbineModelId,
      });
      toast.success('DLC definition created');
      setShowCreateModal(false);
      setNewName('');
      setNewTurbineModelId('');
      const data = await loadData();
      if (data) {
        const found = data.find((d) => d.id === created.id);
        if (found) {
          setSelectedDlc(found);
          setEditForm(dlcToEditable(found));
        }
      }
    } catch (err) {
      toast.error('Failed to create DLC definition');
    } finally {
      setIsCreating(false);
    }
  };

  const handleSave = async () => {
    if (!projectId || !selectedDlc || !editForm) return;
    setIsSaving(true);
    try {
      const payload: Partial<DLCDefinitionCreate> = {
        name: editForm.name,
        turbine_model_id: editForm.turbine_model_id,
        dlc_cases: editForm.dlc_cases,
        turbsim_params: editForm.turbsim_params,
      };
      const updated = await dlcDefinitionsApi.update(projectId, selectedDlc.id, payload);
      toast.success('DLC definition saved');
      setSelectedDlc(updated);
      setEditForm(dlcToEditable(updated));
      await loadData();
    } catch (err) {
      toast.error('Failed to save DLC definition');
    } finally {
      setIsSaving(false);
    }
  };

  const handleDelete = async () => {
    if (!projectId || !selectedDlc) return;
    if (!window.confirm(`Delete DLC definition "${selectedDlc.name}"?`)) return;
    try {
      await dlcDefinitionsApi.delete(projectId, selectedDlc.id);
      toast.success('DLC definition deleted');
      setSelectedDlc(null);
      setEditForm(null);
      await loadData();
    } catch (err) {
      toast.error('Failed to delete DLC definition');
    }
  };

  const updateCase = (index: number, field: keyof DLCCaseSpec, value: any) => {
    if (!editForm) return;
    const cases = [...editForm.dlc_cases];
    cases[index] = { ...cases[index], [field]: value };
    setEditForm({ ...editForm, dlc_cases: cases });
  };

  const addCase = () => {
    if (!editForm) return;
    setEditForm({ ...editForm, dlc_cases: [...editForm.dlc_cases, { ...DEFAULT_CASE }] });
  };

  const removeCase = (index: number) => {
    if (!editForm) return;
    setEditForm({ ...editForm, dlc_cases: editForm.dlc_cases.filter((_, i) => i !== index) });
  };

  const updateTurbsim = (field: keyof TurbSimParams, value: any) => {
    if (!editForm) return;
    setEditForm({
      ...editForm,
      turbsim_params: { ...editForm.turbsim_params, [field]: value },
    });
  };

  const parseNumberArray = (str: string): number[] => {
    return str
      .split(',')
      .map((s) => s.trim())
      .filter((s) => s !== '')
      .map(Number)
      .filter((n) => !isNaN(n));
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
          <h2 className="text-xl font-bold text-slate-800">DLC Matrix</h2>
          <p className="text-sm text-slate-500">
            Define Design Load Case sets for IEC certification
          </p>
        </div>
        <button
          onClick={() => setShowCreateModal(true)}
          className="btn-primary flex items-center gap-2"
        >
          <Plus className="h-4 w-4" />
          New DLC Definition
        </button>
      </div>

      {/* Create Modal */}
      {showCreateModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40">
          <div className="bg-white rounded-xl shadow-xl p-6 w-full max-w-md">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-lg font-semibold text-slate-800">New DLC Definition</h3>
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
                  placeholder="e.g., IEC Class I DLC Set"
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                />
              </div>
              <div>
                <label className="block text-sm font-medium text-slate-700 mb-1">Turbine Model</label>
                <select
                  value={newTurbineModelId}
                  onChange={(e) => setNewTurbineModelId(e.target.value)}
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                >
                  <option value="">-- Select --</option>
                  {turbineModels.map((m) => (
                    <option key={m.id} value={m.id}>{m.name}</option>
                  ))}
                </select>
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
                disabled={isCreating || !newName.trim() || !newTurbineModelId}
                className="px-4 py-2 text-sm font-medium text-white bg-accent-600 rounded-lg hover:bg-accent-700 disabled:opacity-50 flex items-center gap-2"
              >
                {isCreating && <Loader2 className="h-4 w-4 animate-spin" />}
                Create
              </button>
            </div>
          </div>
        </div>
      )}

      {dlcDefs.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Table2 className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">No DLC definitions</h3>
          <p className="mt-1 text-sm text-slate-500 max-w-sm text-center">
            Create Design Load Case definitions to configure the simulation matrix.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* DLC list */}
          <div className="space-y-2">
            {dlcDefs.map((dlc) => (
              <button
                key={dlc.id}
                onClick={() => handleSelect(dlc)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedDlc?.id === dlc.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <span className="font-medium text-slate-800">{dlc.name}</span>
                <div className="mt-1 text-xs text-slate-500">
                  {dlc.total_case_count} total cases &middot; {dlc.status}
                </div>
              </button>
            ))}
          </div>

          {/* DLC detail */}
          {selectedDlc && editForm && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6 space-y-6">
              <div className="flex items-center justify-between">
                <h3 className="text-lg font-semibold text-slate-800">Edit DLC Definition</h3>
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

              {/* Name + turbine model */}
              <div className="grid grid-cols-2 gap-4">
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Name</label>
                  <input
                    type="text"
                    value={editForm.name}
                    onChange={(e) => setEditForm({ ...editForm, name: e.target.value })}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">Turbine Model</label>
                  <select
                    value={editForm.turbine_model_id}
                    onChange={(e) => setEditForm({ ...editForm, turbine_model_id: e.target.value })}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  >
                    <option value="">-- Select --</option>
                    {turbineModels.map((m) => (
                      <option key={m.id} value={m.id}>{m.name}</option>
                    ))}
                  </select>
                </div>
              </div>

              {/* DLC Cases */}
              <div>
                <div className="flex items-center justify-between mb-3">
                  <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider">
                    DLC Cases ({editForm.dlc_cases.length})
                  </h4>
                  <button
                    onClick={addCase}
                    className="flex items-center gap-1 px-3 py-1.5 text-xs font-medium text-accent-700 bg-accent-50 rounded-lg hover:bg-accent-100"
                  >
                    <Plus className="h-3.5 w-3.5" />
                    Add DLC Case
                  </button>
                </div>

                {editForm.dlc_cases.length > 0 ? (
                  <div className="space-y-3">
                    {editForm.dlc_cases.map((c, i) => (
                      <div key={i} className="rounded-lg border border-slate-200 p-4">
                        <div className="flex items-center justify-between mb-3">
                          <span className="text-sm font-medium text-slate-700">Case #{i + 1}</span>
                          <button onClick={() => removeCase(i)} className="text-slate-400 hover:text-danger-500">
                            <Trash2 className="h-3.5 w-3.5" />
                          </button>
                        </div>
                        <div className="grid grid-cols-2 gap-3">
                          <div>
                            <label className="block text-xs text-slate-500 mb-0.5">DLC Number</label>
                            <input
                              type="text"
                              value={c.dlc_number}
                              onChange={(e) => updateCase(i, 'dlc_number', e.target.value)}
                              className="w-full rounded border border-slate-200 px-2 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </div>
                          <div>
                            <label className="block text-xs text-slate-500 mb-0.5">Seeds</label>
                            <input
                              type="number"
                              value={c.seeds}
                              onChange={(e) => updateCase(i, 'seeds', parseInt(e.target.value) || 1)}
                              className="w-full rounded border border-slate-200 px-2 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </div>
                          <div>
                            <label className="block text-xs text-slate-500 mb-0.5">Wind Speeds (comma-separated)</label>
                            <input
                              type="text"
                              value={c.wind_speeds.join(', ')}
                              onChange={(e) => updateCase(i, 'wind_speeds', parseNumberArray(e.target.value))}
                              className="w-full rounded border border-slate-200 px-2 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </div>
                          <div>
                            <label className="block text-xs text-slate-500 mb-0.5">Yaw Misalignments (comma-separated)</label>
                            <input
                              type="text"
                              value={c.yaw_misalignments.join(', ')}
                              onChange={(e) => updateCase(i, 'yaw_misalignments', parseNumberArray(e.target.value))}
                              className="w-full rounded border border-slate-200 px-2 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                            />
                          </div>
                        </div>
                      </div>
                    ))}
                  </div>
                ) : (
                  <div className="rounded-lg border-2 border-dashed border-slate-200 p-8 text-center">
                    <p className="text-sm text-slate-500">No DLC cases. Click "Add DLC Case" to begin.</p>
                  </div>
                )}
              </div>

              {/* Total case count */}
              <div className="rounded-lg bg-slate-50 p-3">
                <p className="text-sm text-slate-600">
                  Total simulation cases:{' '}
                  <span className="font-semibold text-slate-800">
                    {computeTotalCases(editForm.dlc_cases)}
                  </span>
                </p>
              </div>

              {/* TurbSim params */}
              <div>
                <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider mb-3">
                  TurbSim Parameters
                </h4>
                <div className="grid grid-cols-2 gap-4">
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">Turbulence Model</label>
                    <input type="text" value={editForm.turbsim_params.turbulence_model}
                      onChange={(e) => updateTurbsim('turbulence_model', e.target.value)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">IEC Standard</label>
                    <input type="text" value={editForm.turbsim_params.iec_standard}
                      onChange={(e) => updateTurbsim('iec_standard', e.target.value)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">IEC Turbulence Class</label>
                    <input type="text" value={editForm.turbsim_params.iec_turbc}
                      onChange={(e) => updateTurbsim('iec_turbc', e.target.value)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">Grid Height (m)</label>
                    <input type="number" value={editForm.turbsim_params.grid_height}
                      onChange={(e) => updateTurbsim('grid_height', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">Grid Width (m)</label>
                    <input type="number" value={editForm.turbsim_params.grid_width}
                      onChange={(e) => updateTurbsim('grid_width', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">Num Grid Z</label>
                    <input type="number" value={editForm.turbsim_params.num_grid_z}
                      onChange={(e) => updateTurbsim('num_grid_z', parseInt(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">Num Grid Y</label>
                    <input type="number" value={editForm.turbsim_params.num_grid_y}
                      onChange={(e) => updateTurbsim('num_grid_y', parseInt(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">Time Step (s)</label>
                    <input type="number" step="0.01" value={editForm.turbsim_params.time_step}
                      onChange={(e) => updateTurbsim('time_step', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">Analysis Time (s)</label>
                    <input type="number" value={editForm.turbsim_params.analysis_time}
                      onChange={(e) => updateTurbsim('analysis_time', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                  <div>
                    <label className="block text-xs text-slate-500 mb-0.5">Ref Height (m)</label>
                    <input type="number" value={editForm.turbsim_params.ref_height}
                      onChange={(e) => updateTurbsim('ref_height', parseFloat(e.target.value) || 0)}
                      className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500" />
                  </div>
                </div>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
