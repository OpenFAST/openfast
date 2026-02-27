import { useEffect, useState, useCallback, useRef } from 'react';
import { useParams } from 'react-router-dom';
import { simulationsApi, dlcDefinitionsApi, turbineModelsApi } from '@/api/client';
import type { Simulation, SimulationCase, DLCDefinition, TurbineModel } from '@/types';
import type { SimulationCreate } from '@/api/client';
import {
  Play,
  Plus,
  Loader2,
  StopCircle,
  CheckCircle2,
  XCircle,
  Clock,
  X,
} from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';

export default function SimulationRunner() {
  const { projectId } = useParams<{ projectId: string }>();
  const [simulations, setSimulations] = useState<Simulation[]>([]);
  const [dlcDefs, setDlcDefs] = useState<DLCDefinition[]>([]);
  const [turbineModels, setTurbineModels] = useState<TurbineModel[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedSim, setSelectedSim] = useState<Simulation | null>(null);
  const [cases, setCases] = useState<SimulationCase[]>([]);
  const [showCreateModal, setShowCreateModal] = useState(false);
  const [newName, setNewName] = useState('');
  const [newDlcId, setNewDlcId] = useState('');
  const [newTurbineModelId, setNewTurbineModelId] = useState('');
  const [isCreating, setIsCreating] = useState(false);
  const refreshRef = useRef<ReturnType<typeof setInterval> | null>(null);

  const loadSimulations = useCallback(async () => {
    if (!projectId) return;
    try {
      const data = await simulationsApi.list(projectId);
      setSimulations(data);
      return data;
    } catch (err) {
      toast.error('Failed to load simulations');
      return [];
    }
  }, [projectId]);

  const loadSupportData = useCallback(async () => {
    if (!projectId) return;
    try {
      const [defs, models] = await Promise.all([
        dlcDefinitionsApi.list(projectId),
        turbineModelsApi.list(projectId),
      ]);
      setDlcDefs(defs);
      setTurbineModels(models);
    } catch {
      // non-critical
    }
  }, [projectId]);

  useEffect(() => {
    setIsLoading(true);
    Promise.all([loadSimulations(), loadSupportData()])
      .then(([sims]) => {
        if (sims && sims.length > 0) {
          setSelectedSim(sims[0]);
        }
      })
      .finally(() => setIsLoading(false));
  }, [loadSimulations, loadSupportData]);

  // Load cases when selection changes
  useEffect(() => {
    if (!projectId || !selectedSim) {
      setCases([]);
      return;
    }
    simulationsApi
      .getCases(projectId, selectedSim.id)
      .then(setCases)
      .catch(() => setCases([]));
  }, [projectId, selectedSim?.id]);

  // Auto-refresh while running
  useEffect(() => {
    if (refreshRef.current) {
      clearInterval(refreshRef.current);
      refreshRef.current = null;
    }

    if (selectedSim?.status === 'running' && projectId) {
      refreshRef.current = setInterval(async () => {
        try {
          const [updatedSim, updatedCases] = await Promise.all([
            simulationsApi.get(projectId, selectedSim.id),
            simulationsApi.getCases(projectId, selectedSim.id),
          ]);
          setSelectedSim(updatedSim);
          setCases(updatedCases);
          setSimulations((prev) =>
            prev.map((s) => (s.id === updatedSim.id ? updatedSim : s)),
          );
          // Stop refreshing if no longer running
          if (updatedSim.status !== 'running' && refreshRef.current) {
            clearInterval(refreshRef.current);
            refreshRef.current = null;
          }
        } catch {
          // ignore refresh errors
        }
      }, 3000);
    }

    return () => {
      if (refreshRef.current) {
        clearInterval(refreshRef.current);
        refreshRef.current = null;
      }
    };
  }, [selectedSim?.id, selectedSim?.status, projectId]);

  const handleSelect = (sim: Simulation) => {
    setSelectedSim(sim);
  };

  const handleCreate = async () => {
    if (!projectId || !newName.trim() || !newDlcId || !newTurbineModelId) return;
    setIsCreating(true);
    try {
      const payload: SimulationCreate = {
        name: newName.trim(),
        dlc_definition_id: newDlcId,
        turbine_model_id: newTurbineModelId,
      };
      const created = await simulationsApi.create(projectId, payload);
      toast.success('Simulation created');
      setShowCreateModal(false);
      setNewName('');
      setNewDlcId('');
      setNewTurbineModelId('');
      const data = await loadSimulations();
      if (data) {
        const found = data.find((s) => s.id === created.id);
        if (found) setSelectedSim(found);
      }
    } catch (err) {
      toast.error('Failed to create simulation');
    } finally {
      setIsCreating(false);
    }
  };

  const handleStart = async () => {
    if (!projectId || !selectedSim) return;
    try {
      const updated = await simulationsApi.start(projectId, selectedSim.id);
      toast.success('Simulation started');
      setSelectedSim(updated);
      setSimulations((prev) =>
        prev.map((s) => (s.id === updated.id ? updated : s)),
      );
    } catch (err) {
      toast.error('Failed to start simulation');
    }
  };

  const handleCancel = async () => {
    if (!projectId || !selectedSim) return;
    try {
      const updated = await simulationsApi.cancel(projectId, selectedSim.id);
      toast.success('Simulation cancelled');
      setSelectedSim(updated);
      setSimulations((prev) =>
        prev.map((s) => (s.id === updated.id ? updated : s)),
      );
    } catch (err) {
      toast.error('Failed to cancel simulation');
    }
  };

  const statusIcon = (status: string) => {
    switch (status) {
      case 'completed':
        return <CheckCircle2 className="h-4 w-4 text-success-500" />;
      case 'failed':
        return <XCircle className="h-4 w-4 text-danger-500" />;
      case 'running':
        return <Loader2 className="h-4 w-4 text-accent-500 animate-spin" />;
      case 'cancelled':
        return <StopCircle className="h-4 w-4 text-warning-500" />;
      default:
        return <Clock className="h-4 w-4 text-slate-400" />;
    }
  };

  const statusBadge = (status: string) => {
    const colors: Record<string, string> = {
      pending: 'bg-slate-100 text-slate-700',
      running: 'bg-accent-100 text-accent-700',
      completed: 'bg-success-100 text-success-700',
      failed: 'bg-danger-100 text-danger-700',
      cancelled: 'bg-warning-100 text-warning-700',
    };
    return (
      <span className={clsx('inline-flex items-center gap-1 rounded-full px-2.5 py-0.5 text-xs font-medium', colors[status] || colors.pending)}>
        {statusIcon(status)}
        {status}
      </span>
    );
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
          <h2 className="text-xl font-bold text-slate-800">Simulation Runner</h2>
          <p className="text-sm text-slate-500">
            Launch and monitor OpenFAST simulations
          </p>
        </div>
        <button
          onClick={() => setShowCreateModal(true)}
          className="btn-primary flex items-center gap-2"
        >
          <Plus className="h-4 w-4" />
          New Simulation
        </button>
      </div>

      {/* Create Modal */}
      {showCreateModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/40">
          <div className="bg-white rounded-xl shadow-xl p-6 w-full max-w-md">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-lg font-semibold text-slate-800">New Simulation</h3>
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
                  placeholder="e.g., Run 001"
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                />
              </div>
              <div>
                <label className="block text-sm font-medium text-slate-700 mb-1">DLC Definition</label>
                <select
                  value={newDlcId}
                  onChange={(e) => setNewDlcId(e.target.value)}
                  className="w-full rounded-lg border border-slate-300 px-3 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
                >
                  <option value="">-- Select --</option>
                  {dlcDefs.map((d) => (
                    <option key={d.id} value={d.id}>{d.name}</option>
                  ))}
                </select>
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
                disabled={isCreating || !newName.trim() || !newDlcId || !newTurbineModelId}
                className="px-4 py-2 text-sm font-medium text-white bg-accent-600 rounded-lg hover:bg-accent-700 disabled:opacity-50 flex items-center gap-2"
              >
                {isCreating && <Loader2 className="h-4 w-4 animate-spin" />}
                Create
              </button>
            </div>
          </div>
        </div>
      )}

      {simulations.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Play className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">No simulations</h3>
          <p className="mt-1 text-sm text-slate-500 max-w-sm text-center">
            Create a simulation by selecting a turbine model and DLC definition.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Simulation list */}
          <div className="space-y-2">
            {simulations.map((sim) => (
              <button
                key={sim.id}
                onClick={() => handleSelect(sim)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedSim?.id === sim.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <div className="flex items-center justify-between">
                  <span className="font-medium text-slate-800">{sim.name}</span>
                  {statusIcon(sim.status)}
                </div>
                <div className="mt-1 text-xs text-slate-500">
                  {sim.completed_cases}/{sim.total_cases} cases &middot; {sim.status}
                </div>
                {sim.status === 'running' && (
                  <div className="mt-2 h-1.5 w-full rounded-full bg-slate-200">
                    <div
                      className="h-full rounded-full bg-accent-500 transition-all duration-500"
                      style={{
                        width: `${sim.total_cases > 0 ? (sim.completed_cases / sim.total_cases) * 100 : 0}%`,
                      }}
                    />
                  </div>
                )}
              </button>
            ))}
          </div>

          {/* Simulation detail */}
          {selectedSim && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-6">
                <div>
                  <h3 className="text-lg font-semibold text-slate-800">{selectedSim.name}</h3>
                  <div className="mt-1">{statusBadge(selectedSim.status)}</div>
                </div>
                <div className="flex items-center gap-2">
                  {selectedSim.status === 'pending' && (
                    <button
                      onClick={handleStart}
                      className="flex items-center gap-2 px-4 py-2 text-sm font-medium text-white bg-accent-600 rounded-lg hover:bg-accent-700"
                    >
                      <Play className="h-4 w-4" />
                      Start
                    </button>
                  )}
                  {selectedSim.status === 'running' && (
                    <button
                      onClick={handleCancel}
                      className="flex items-center gap-2 px-4 py-2 text-sm font-medium text-white bg-danger-600 rounded-lg hover:bg-danger-700"
                    >
                      <StopCircle className="h-4 w-4" />
                      Cancel
                    </button>
                  )}
                </div>
              </div>

              {/* Overall progress */}
              <div className="mb-6 rounded-lg bg-slate-50 p-4">
                <div className="flex items-center justify-between mb-2">
                  <span className="text-sm font-medium text-slate-700">Overall Progress</span>
                  <span className="text-sm font-mono text-slate-600">
                    {selectedSim.completed_cases}/{selectedSim.total_cases}
                  </span>
                </div>
                <div className="h-2 w-full rounded-full bg-slate-200">
                  <div
                    className="h-full rounded-full bg-accent-500 transition-all duration-500"
                    style={{
                      width: `${selectedSim.total_cases > 0 ? (selectedSim.completed_cases / selectedSim.total_cases) * 100 : 0}%`,
                    }}
                  />
                </div>
                {selectedSim.failed_cases > 0 && (
                  <p className="mt-2 text-xs text-danger-600">
                    {selectedSim.failed_cases} case{selectedSim.failed_cases !== 1 ? 's' : ''} failed
                  </p>
                )}
              </div>

              {/* Cases table */}
              {cases.length > 0 && (
                <div className="overflow-x-auto rounded-lg border border-slate-200">
                  <table className="min-w-full text-sm">
                    <thead className="bg-slate-50">
                      <tr>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">DLC</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Wind (m/s)</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Seed</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Status</th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Progress</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-100">
                      {cases.map((c) => (
                        <tr key={c.id} className="hover:bg-slate-50">
                          <td className="px-3 py-2 font-medium text-slate-800">{c.dlc_number}</td>
                          <td className="px-3 py-2 text-slate-700">{c.wind_speed}</td>
                          <td className="px-3 py-2 text-slate-700">{c.seed_number}</td>
                          <td className="px-3 py-2">
                            <div className="flex items-center gap-1.5">
                              {statusIcon(c.status)}
                              <span className="text-slate-700 capitalize">{c.status}</span>
                            </div>
                          </td>
                          <td className="px-3 py-2">
                            <div className="flex items-center gap-2">
                              <div className="h-1.5 w-20 rounded-full bg-slate-200">
                                <div
                                  className={clsx(
                                    'h-full rounded-full transition-all duration-300',
                                    c.status === 'completed'
                                      ? 'bg-success-500'
                                      : c.status === 'failed'
                                        ? 'bg-danger-500'
                                        : 'bg-accent-500',
                                  )}
                                  style={{ width: `${c.progress_percent}%` }}
                                />
                              </div>
                              <span className="text-xs font-mono text-slate-500 w-8">
                                {c.progress_percent}%
                              </span>
                            </div>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
