import { useEffect, useState, useCallback } from 'react';
import { useParams } from 'react-router-dom';
import { controllersApi } from '@/api/client';
import type { Controller, ControllerCreate } from '@/types';
import { Gauge, Plus, Loader2, Trash2, Save, X } from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';

interface EditableController {
  name: string;
  controller_type: string;
  pcmode: number;
  vscontrl: number;
  parameters: Record<string, any>;
  dll_filename: string;
  dll_procname: string;
}

function controllerToEditable(ctrl: Controller): EditableController {
  return {
    name: ctrl.name,
    controller_type: ctrl.controller_type,
    pcmode: ctrl.pcmode,
    vscontrl: ctrl.vscontrl,
    parameters: ctrl.parameters ? { ...ctrl.parameters } : {
      blade_pitch_min: 0,
      blade_pitch_max: 90,
      pitch_rate_limit: 8,
      kp_pitch: 0,
      ki_pitch: 0,
      rated_gen_speed: 0,
      rgn2k: 0,
      rated_torque: 0,
      min_gen_torque: 0,
      max_gen_torque: 0,
      gen_model: 1,
      gen_eff: 0.944,
    },
    dll_filename: ctrl.dll_filename || '',
    dll_procname: ctrl.dll_procname || '',
  };
}

const PITCH_PARAMS = ['blade_pitch_min', 'blade_pitch_max', 'pitch_rate_limit', 'kp_pitch', 'ki_pitch'];
const TORQUE_PARAMS = ['rated_gen_speed', 'rgn2k', 'rated_torque', 'min_gen_torque', 'max_gen_torque'];
const GENERAL_PARAMS = ['gen_model', 'gen_eff'];

const PARAM_LABELS: Record<string, string> = {
  blade_pitch_min: 'Blade Pitch Min (deg)',
  blade_pitch_max: 'Blade Pitch Max (deg)',
  pitch_rate_limit: 'Pitch Rate Limit (deg/s)',
  kp_pitch: 'Kp Pitch',
  ki_pitch: 'Ki Pitch',
  rated_gen_speed: 'Rated Gen Speed (rpm)',
  rgn2k: 'Region 2 K',
  rated_torque: 'Rated Torque (Nm)',
  min_gen_torque: 'Min Gen Torque (Nm)',
  max_gen_torque: 'Max Gen Torque (Nm)',
  gen_model: 'Generator Model',
  gen_eff: 'Generator Efficiency',
};

export default function ControllerDesigner() {
  const { projectId } = useParams<{ projectId: string }>();
  const [controllers, setControllers] = useState<Controller[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedController, setSelectedController] = useState<Controller | null>(null);
  const [editForm, setEditForm] = useState<EditableController | null>(null);
  const [isSaving, setIsSaving] = useState(false);
  const [showCreateModal, setShowCreateModal] = useState(false);
  const [newName, setNewName] = useState('');
  const [newType, setNewType] = useState('baseline');
  const [isCreating, setIsCreating] = useState(false);

  const loadControllers = useCallback(async () => {
    if (!projectId) return;
    try {
      const data = await controllersApi.list(projectId);
      setControllers(data);
      return data;
    } catch (err) {
      toast.error('Failed to load controllers');
      return [];
    }
  }, [projectId]);

  useEffect(() => {
    setIsLoading(true);
    loadControllers()
      .then((data) => {
        if (data && data.length > 0) {
          setSelectedController(data[0]);
          setEditForm(controllerToEditable(data[0]));
        }
      })
      .finally(() => setIsLoading(false));
  }, [loadControllers]);

  const handleSelect = (ctrl: Controller) => {
    setSelectedController(ctrl);
    setEditForm(controllerToEditable(ctrl));
  };

  const handleCreate = async () => {
    if (!projectId || !newName.trim()) return;
    setIsCreating(true);
    try {
      const created = await controllersApi.create(projectId, {
        name: newName.trim(),
        controller_type: newType,
      } as ControllerCreate);
      toast.success('Controller created');
      setShowCreateModal(false);
      setNewName('');
      setNewType('baseline');
      const data = await loadControllers();
      if (data) {
        const found = data.find((c) => c.id === created.id);
        if (found) {
          setSelectedController(found);
          setEditForm(controllerToEditable(found));
        }
      }
    } catch (err) {
      toast.error('Failed to create controller');
    } finally {
      setIsCreating(false);
    }
  };

  const handleSave = async () => {
    if (!projectId || !selectedController || !editForm) return;
    setIsSaving(true);
    try {
      const payload: Partial<ControllerCreate> = {
        name: editForm.name,
        controller_type: editForm.controller_type,
        pcmode: editForm.pcmode,
        vscontrl: editForm.vscontrl,
        parameters: editForm.parameters,
        dll_filename: editForm.dll_filename || undefined,
        dll_procname: editForm.dll_procname || undefined,
      };
      const updated = await controllersApi.update(projectId, selectedController.id, payload);
      toast.success('Controller saved');
      setSelectedController(updated);
      setEditForm(controllerToEditable(updated));
      await loadControllers();
    } catch (err) {
      toast.error('Failed to save controller');
    } finally {
      setIsSaving(false);
    }
  };

  const handleDelete = async () => {
    if (!projectId || !selectedController) return;
    if (!window.confirm(`Delete controller "${selectedController.name}"?`)) return;
    try {
      await controllersApi.delete(projectId, selectedController.id);
      toast.success('Controller deleted');
      setSelectedController(null);
      setEditForm(null);
      await loadControllers();
    } catch (err) {
      toast.error('Failed to delete controller');
    }
  };

  const updateField = (field: keyof EditableController, value: any) => {
    if (!editForm) return;
    setEditForm({ ...editForm, [field]: value });
  };

  const updateParam = (key: string, value: number) => {
    if (!editForm) return;
    setEditForm({
      ...editForm,
      parameters: { ...editForm.parameters, [key]: value },
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
          <h2 className="text-xl font-bold text-slate-800">Controller Designer</h2>
          <p className="text-sm text-slate-500">
            Configure pitch and torque controller parameters
          </p>
        </div>
        <button
          onClick={() => setShowCreateModal(true)}
          className="btn-primary flex items-center gap-2"
        >
          <Plus className="h-4 w-4" />
          New Controller
        </button>
      </div>

      {/* Create Modal */}
      {showCreateModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40">
          <div className="bg-white rounded-xl shadow-xl p-6 w-full max-w-md">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-lg font-semibold text-slate-800">New Controller</h3>
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
                  placeholder="e.g., Baseline Controller"
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                />
              </div>
              <div>
                <label className="block text-sm font-medium text-slate-700 mb-1">Controller Type</label>
                <select
                  value={newType}
                  onChange={(e) => setNewType(e.target.value)}
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                >
                  <option value="baseline">Baseline</option>
                  <option value="ROSCO">ROSCO</option>
                  <option value="custom">Custom</option>
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

      {controllers.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Gauge className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">No controllers defined</h3>
          <p className="mt-1 text-sm text-slate-500">
            Create your first controller configuration.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Controller list */}
          <div className="space-y-2">
            {controllers.map((ctrl) => (
              <button
                key={ctrl.id}
                onClick={() => handleSelect(ctrl)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedController?.id === ctrl.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <span className="font-medium text-slate-800">{ctrl.name}</span>
                <div className="mt-1 text-xs text-slate-500">
                  {ctrl.controller_type} &middot; PCMode {ctrl.pcmode}
                </div>
              </button>
            ))}
          </div>

          {/* Controller detail */}
          {selectedController && editForm && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-6">
                <h3 className="text-lg font-semibold text-slate-800">Edit Controller</h3>
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
                  <label className="block text-xs font-medium text-slate-500 mb-1">Controller Type</label>
                  <select
                    value={editForm.controller_type}
                    onChange={(e) => updateField('controller_type', e.target.value)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  >
                    <option value="baseline">Baseline</option>
                    <option value="ROSCO">ROSCO</option>
                    <option value="custom">Custom</option>
                  </select>
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">PCMode</label>
                  <input
                    type="number"
                    value={editForm.pcmode}
                    onChange={(e) => updateField('pcmode', parseInt(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">VSContrl</label>
                  <input
                    type="number"
                    value={editForm.vscontrl}
                    onChange={(e) => updateField('vscontrl', parseInt(e.target.value) || 0)}
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
              </div>

              {/* Parameters */}
              <div className="grid grid-cols-2 gap-6 mb-6">
                {/* Pitch Parameters */}
                <div className="space-y-3">
                  <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider">
                    Pitch Control
                  </h4>
                  <div className="space-y-2">
                    {PITCH_PARAMS.map((key) => (
                      <div key={key}>
                        <label className="block text-xs text-slate-500 mb-0.5">{PARAM_LABELS[key]}</label>
                        <input
                          type="number"
                          step="any"
                          value={editForm.parameters[key] ?? 0}
                          onChange={(e) => updateParam(key, parseFloat(e.target.value) || 0)}
                          className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm font-mono focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                        />
                      </div>
                    ))}
                  </div>
                </div>

                {/* Torque Parameters */}
                <div className="space-y-3">
                  <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider">
                    Torque Control
                  </h4>
                  <div className="space-y-2">
                    {TORQUE_PARAMS.map((key) => (
                      <div key={key}>
                        <label className="block text-xs text-slate-500 mb-0.5">{PARAM_LABELS[key]}</label>
                        <input
                          type="number"
                          step="any"
                          value={editForm.parameters[key] ?? 0}
                          onChange={(e) => updateParam(key, parseFloat(e.target.value) || 0)}
                          className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm font-mono focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                        />
                      </div>
                    ))}
                  </div>
                </div>
              </div>

              {/* General Parameters */}
              <div className="mb-6">
                <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider mb-3">
                  General
                </h4>
                <div className="grid grid-cols-2 gap-4">
                  {GENERAL_PARAMS.map((key) => (
                    <div key={key}>
                      <label className="block text-xs text-slate-500 mb-0.5">{PARAM_LABELS[key]}</label>
                      <input
                        type="number"
                        step="any"
                        value={editForm.parameters[key] ?? 0}
                        onChange={(e) => updateParam(key, parseFloat(e.target.value) || 0)}
                        className="w-full rounded-lg border border-slate-300 px-3 py-1.5 text-sm font-mono focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                      />
                    </div>
                  ))}
                </div>
              </div>

              {/* DLL fields */}
              <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider mb-3">
                DLL Configuration
              </h4>
              <div className="grid grid-cols-2 gap-4">
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">DLL Filename</label>
                  <input
                    type="text"
                    value={editForm.dll_filename}
                    onChange={(e) => updateField('dll_filename', e.target.value)}
                    placeholder="e.g., DISCON.dll"
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
                <div>
                  <label className="block text-xs font-medium text-slate-500 mb-1">DLL Procedure Name</label>
                  <input
                    type="text"
                    value={editForm.dll_procname}
                    onChange={(e) => updateField('dll_procname', e.target.value)}
                    placeholder="e.g., DISCON"
                    className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                  />
                </div>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
