import { useEffect, useState, type FormEvent } from 'react';
import { useNavigate } from 'react-router-dom';
import { useProjectStore } from '@/stores/projectStore';
import type { ProjectCreate, WindClass, TurbulenceClass } from '@/types';
import {
  Plus,
  Wind,
  Calendar,
  Ruler,
  Zap,
  X,
  Loader2,
  FolderOpen,
} from 'lucide-react';
import clsx from 'clsx';

const WIND_CLASSES: WindClass[] = ['I', 'II', 'III', 'S'];
const TURBULENCE_CLASSES: TurbulenceClass[] = ['A', 'B', 'C'];

const defaultFormData: ProjectCreate = {
  name: '',
  description: '',
  wind_class: 'II',
  turbulence_class: 'B',
  rated_power_kw: 5000,
  rotor_diameter_m: 126,
  hub_height_m: 90,
  cut_in_speed: 3,
  rated_speed: 11.4,
  cut_out_speed: 25,
  dt: 0.0125,
  t_max: 660,
};

export default function Dashboard() {
  const projects = useProjectStore((s) => s.projects);
  const isLoading = useProjectStore((s) => s.isLoading);
  const fetchProjects = useProjectStore((s) => s.fetchProjects);
  const createProject = useProjectStore((s) => s.createProject);
  const navigate = useNavigate();

  const [showModal, setShowModal] = useState(false);
  const [formData, setFormData] = useState<ProjectCreate>(defaultFormData);
  const [isCreating, setIsCreating] = useState(false);

  useEffect(() => {
    fetchProjects();
  }, [fetchProjects]);

  const handleCreateProject = async (e: FormEvent) => {
    e.preventDefault();
    if (!formData.name.trim()) return;

    setIsCreating(true);
    try {
      const project = await createProject(formData);
      setShowModal(false);
      setFormData(defaultFormData);
      navigate(`/projects/${project.id}/tower`);
    } catch {
      // Error handled in store
    } finally {
      setIsCreating(false);
    }
  };

  const updateForm = <K extends keyof ProjectCreate>(
    key: K,
    value: ProjectCreate[K],
  ) => {
    setFormData((prev) => ({ ...prev, [key]: value }));
  };

  return (
    <div className="p-8">
      {/* Header */}
      <div className="mb-8 flex items-center justify-between">
        <div>
          <h1 className="text-2xl font-bold text-slate-800">Your Projects</h1>
          <p className="mt-1 text-sm text-slate-500">
            Manage and design your wind turbine configurations
          </p>
        </div>
        <button onClick={() => setShowModal(true)} className="btn-primary">
          <Plus className="h-4 w-4" />
          New Project
        </button>
      </div>

      {/* Loading state */}
      {isLoading && projects.length === 0 && (
        <div className="flex items-center justify-center py-24">
          <Loader2 className="h-8 w-8 animate-spin text-accent-500" />
        </div>
      )}

      {/* Empty state */}
      {!isLoading && projects.length === 0 && (
        <div className="flex flex-col items-center justify-center rounded-2xl border-2 border-dashed border-slate-200 py-24">
          <div className="mb-4 flex h-16 w-16 items-center justify-center rounded-2xl bg-slate-100">
            <FolderOpen className="h-8 w-8 text-slate-400" />
          </div>
          <h3 className="text-lg font-semibold text-slate-700">
            No projects yet
          </h3>
          <p className="mt-1 text-sm text-slate-500 max-w-sm text-center">
            Create your first wind turbine project to start designing towers,
            blades, and controllers.
          </p>
          <button
            onClick={() => setShowModal(true)}
            className="btn-primary mt-6"
          >
            <Plus className="h-4 w-4" />
            Create your first project
          </button>
        </div>
      )}

      {/* Project grid */}
      {projects.length > 0 && (
        <div className="grid grid-cols-1 gap-5 sm:grid-cols-2 lg:grid-cols-3">
          {projects.map((project) => (
            <button
              key={project.id}
              onClick={() => navigate(`/projects/${project.id}/tower`)}
              className="group rounded-xl border border-slate-200 bg-white p-6 text-left shadow-sm transition-all duration-200 hover:border-accent-300 hover:shadow-md hover:shadow-accent-500/5"
            >
              {/* Project header */}
              <div className="mb-4 flex items-start justify-between">
                <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-accent-50 transition-colors group-hover:bg-accent-100">
                  <Wind className="h-5 w-5 text-accent-600" />
                </div>
                <span
                  className={clsx(
                    'inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-medium',
                    project.status === 'active'
                      ? 'bg-success-50 text-success-700'
                      : project.status === 'completed'
                        ? 'bg-accent-50 text-accent-700'
                        : 'bg-slate-100 text-slate-600',
                  )}
                >
                  {project.status || 'Draft'}
                </span>
              </div>

              {/* Name & description */}
              <h3 className="text-base font-semibold text-slate-800 group-hover:text-accent-700 transition-colors">
                {project.name}
              </h3>
              {project.description && (
                <p className="mt-1 text-sm text-slate-500 line-clamp-2">
                  {project.description}
                </p>
              )}

              {/* Specs */}
              <div className="mt-4 grid grid-cols-2 gap-3">
                <div className="flex items-center gap-1.5 text-xs text-slate-500">
                  <Wind className="h-3.5 w-3.5 text-slate-400" />
                  <span>
                    Class {project.wind_class}
                    {project.turbulence_class}
                  </span>
                </div>
                <div className="flex items-center gap-1.5 text-xs text-slate-500">
                  <Ruler className="h-3.5 w-3.5 text-slate-400" />
                  <span>{project.rotor_diameter_m}m rotor</span>
                </div>
                <div className="flex items-center gap-1.5 text-xs text-slate-500">
                  <Zap className="h-3.5 w-3.5 text-slate-400" />
                  <span>{project.rated_power_kw} kW</span>
                </div>
                <div className="flex items-center gap-1.5 text-xs text-slate-500">
                  <Ruler className="h-3.5 w-3.5 text-slate-400" />
                  <span>{project.hub_height_m}m hub</span>
                </div>
              </div>

              {/* Date */}
              <div className="mt-4 flex items-center gap-1.5 text-xs text-slate-400 border-t border-slate-100 pt-3">
                <Calendar className="h-3.5 w-3.5" />
                <span>
                  Created{' '}
                  {new Date(project.created_at).toLocaleDateString('en-US', {
                    month: 'short',
                    day: 'numeric',
                    year: 'numeric',
                  })}
                </span>
              </div>
            </button>
          ))}
        </div>
      )}

      {/* New Project Modal */}
      {showModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center">
          {/* Backdrop */}
          <div
            className="absolute inset-0 bg-black/50 backdrop-blur-sm"
            onClick={() => setShowModal(false)}
          />

          {/* Modal */}
          <div className="relative z-10 w-full max-w-2xl max-h-[90vh] overflow-y-auto rounded-2xl border border-slate-700/50 bg-surface-dark-secondary p-8 shadow-2xl mx-4">
            {/* Close button */}
            <button
              onClick={() => setShowModal(false)}
              className="absolute right-4 top-4 flex h-8 w-8 items-center justify-center rounded-lg text-slate-400 hover:bg-slate-700 hover:text-slate-200 transition-colors"
            >
              <X className="h-5 w-5" />
            </button>

            <h2 className="text-xl font-bold text-white mb-1">
              Create New Project
            </h2>
            <p className="text-sm text-slate-400 mb-6">
              Define your wind turbine project parameters
            </p>

            <form onSubmit={handleCreateProject} className="space-y-6">
              {/* Basic Info */}
              <div className="space-y-4">
                <h3 className="text-sm font-semibold text-slate-300 uppercase tracking-wider">
                  General
                </h3>

                <div className="grid grid-cols-1 gap-4 sm:grid-cols-2">
                  <div className="sm:col-span-2">
                    <label htmlFor="proj-name" className="label">
                      Project name
                    </label>
                    <input
                      id="proj-name"
                      type="text"
                      value={formData.name}
                      onChange={(e) => updateForm('name', e.target.value)}
                      className="input-field"
                      placeholder="e.g., NREL 5MW Reference"
                      required
                      autoFocus
                    />
                  </div>

                  <div className="sm:col-span-2">
                    <label htmlFor="proj-desc" className="label">
                      Description{' '}
                      <span className="text-slate-500 font-normal">
                        (optional)
                      </span>
                    </label>
                    <textarea
                      id="proj-desc"
                      value={formData.description || ''}
                      onChange={(e) =>
                        updateForm('description', e.target.value)
                      }
                      className="input-field resize-none"
                      rows={2}
                      placeholder="Brief project description..."
                    />
                  </div>

                  <div>
                    <label htmlFor="wind-class" className="label">
                      IEC Wind Class
                    </label>
                    <select
                      id="wind-class"
                      value={formData.wind_class}
                      onChange={(e) =>
                        updateForm('wind_class', e.target.value as WindClass)
                      }
                      className="input-field"
                    >
                      {WIND_CLASSES.map((wc) => (
                        <option key={wc} value={wc}>
                          Class {wc}
                        </option>
                      ))}
                    </select>
                  </div>

                  <div>
                    <label htmlFor="turb-class" className="label">
                      Turbulence Class
                    </label>
                    <select
                      id="turb-class"
                      value={formData.turbulence_class}
                      onChange={(e) =>
                        updateForm(
                          'turbulence_class',
                          e.target.value as TurbulenceClass,
                        )
                      }
                      className="input-field"
                    >
                      {TURBULENCE_CLASSES.map((tc) => (
                        <option key={tc} value={tc}>
                          {tc}
                        </option>
                      ))}
                    </select>
                  </div>
                </div>
              </div>

              {/* Turbine parameters */}
              <div className="space-y-4">
                <h3 className="text-sm font-semibold text-slate-300 uppercase tracking-wider">
                  Turbine Parameters
                </h3>

                <div className="grid grid-cols-1 gap-4 sm:grid-cols-3">
                  <div>
                    <label htmlFor="rated-power" className="label">
                      Rated power (kW)
                    </label>
                    <input
                      id="rated-power"
                      type="number"
                      value={formData.rated_power_kw}
                      onChange={(e) =>
                        updateForm('rated_power_kw', Number(e.target.value))
                      }
                      className="input-field"
                      min={0}
                      step={100}
                      required
                    />
                  </div>

                  <div>
                    <label htmlFor="rotor-diam" className="label">
                      Rotor diameter (m)
                    </label>
                    <input
                      id="rotor-diam"
                      type="number"
                      value={formData.rotor_diameter_m}
                      onChange={(e) =>
                        updateForm('rotor_diameter_m', Number(e.target.value))
                      }
                      className="input-field"
                      min={0}
                      step={1}
                      required
                    />
                  </div>

                  <div>
                    <label htmlFor="hub-height" className="label">
                      Hub height (m)
                    </label>
                    <input
                      id="hub-height"
                      type="number"
                      value={formData.hub_height_m}
                      onChange={(e) =>
                        updateForm('hub_height_m', Number(e.target.value))
                      }
                      className="input-field"
                      min={0}
                      step={1}
                      required
                    />
                  </div>
                </div>
              </div>

              {/* Wind speed parameters */}
              <div className="space-y-4">
                <h3 className="text-sm font-semibold text-slate-300 uppercase tracking-wider">
                  Operating Conditions
                </h3>

                <div className="grid grid-cols-1 gap-4 sm:grid-cols-3">
                  <div>
                    <label htmlFor="cut-in" className="label">
                      Cut-in speed (m/s)
                    </label>
                    <input
                      id="cut-in"
                      type="number"
                      value={formData.cut_in_speed}
                      onChange={(e) =>
                        updateForm('cut_in_speed', Number(e.target.value))
                      }
                      className="input-field"
                      min={0}
                      step={0.1}
                      required
                    />
                  </div>

                  <div>
                    <label htmlFor="rated-speed" className="label">
                      Rated speed (m/s)
                    </label>
                    <input
                      id="rated-speed"
                      type="number"
                      value={formData.rated_speed}
                      onChange={(e) =>
                        updateForm('rated_speed', Number(e.target.value))
                      }
                      className="input-field"
                      min={0}
                      step={0.1}
                      required
                    />
                  </div>

                  <div>
                    <label htmlFor="cut-out" className="label">
                      Cut-out speed (m/s)
                    </label>
                    <input
                      id="cut-out"
                      type="number"
                      value={formData.cut_out_speed}
                      onChange={(e) =>
                        updateForm('cut_out_speed', Number(e.target.value))
                      }
                      className="input-field"
                      min={0}
                      step={0.1}
                      required
                    />
                  </div>
                </div>
              </div>

              {/* Simulation defaults */}
              <div className="space-y-4">
                <h3 className="text-sm font-semibold text-slate-300 uppercase tracking-wider">
                  Simulation Defaults
                </h3>

                <div className="grid grid-cols-1 gap-4 sm:grid-cols-2">
                  <div>
                    <label htmlFor="dt" className="label">
                      Time step - dt (s)
                    </label>
                    <input
                      id="dt"
                      type="number"
                      value={formData.dt}
                      onChange={(e) =>
                        updateForm('dt', Number(e.target.value))
                      }
                      className="input-field"
                      min={0.001}
                      step={0.0001}
                      required
                    />
                  </div>

                  <div>
                    <label htmlFor="t-max" className="label">
                      Max time - t_max (s)
                    </label>
                    <input
                      id="t-max"
                      type="number"
                      value={formData.t_max}
                      onChange={(e) =>
                        updateForm('t_max', Number(e.target.value))
                      }
                      className="input-field"
                      min={0}
                      step={10}
                      required
                    />
                  </div>
                </div>
              </div>

              {/* Actions */}
              <div className="flex items-center justify-end gap-3 pt-2 border-t border-slate-700">
                <button
                  type="button"
                  onClick={() => setShowModal(false)}
                  className="btn-secondary"
                >
                  Cancel
                </button>
                <button
                  type="submit"
                  disabled={isCreating || !formData.name.trim()}
                  className="btn-primary"
                >
                  {isCreating ? (
                    <>
                      <Loader2 className="h-4 w-4 animate-spin" />
                      Creating...
                    </>
                  ) : (
                    <>
                      <Plus className="h-4 w-4" />
                      Create Project
                    </>
                  )}
                </button>
              </div>
            </form>
          </div>
        </div>
      )}
    </div>
  );
}
