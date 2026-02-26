import { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { controllersApi } from '@/api/client';
import type { Controller } from '@/types';
import { Gauge, Plus, Loader2, Trash2 } from 'lucide-react';
import clsx from 'clsx';

export default function ControllerDesigner() {
  const { projectId } = useParams<{ projectId: string }>();
  const [controllers, setControllers] = useState<Controller[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedController, setSelectedController] =
    useState<Controller | null>(null);

  useEffect(() => {
    if (!projectId) return;
    setIsLoading(true);
    controllersApi
      .list(projectId)
      .then((data) => {
        setControllers(data);
        if (data.length > 0) setSelectedController(data[0]);
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
          <h2 className="text-xl font-bold text-slate-800">
            Controller Designer
          </h2>
          <p className="text-sm text-slate-500">
            Configure pitch and torque controller parameters
          </p>
        </div>
        <button className="btn-primary">
          <Plus className="h-4 w-4" />
          New Controller
        </button>
      </div>

      {controllers.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Gauge className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">
            No controllers defined
          </h3>
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
                onClick={() => setSelectedController(ctrl)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedController?.id === ctrl.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <span className="font-medium text-slate-800">
                  {ctrl.name}
                </span>
                <div className="mt-1 text-xs text-slate-500">
                  {ctrl.control_mode} &middot; {ctrl.rated_power_kw} kW
                </div>
              </button>
            ))}
          </div>

          {/* Controller detail */}
          {selectedController && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-6">
                <h3 className="text-lg font-semibold text-slate-800">
                  {selectedController.name}
                </h3>
                <button className="text-slate-400 hover:text-danger-500 transition-colors">
                  <Trash2 className="h-4 w-4" />
                </button>
              </div>

              <div className="grid grid-cols-2 gap-6">
                {/* General */}
                <div className="space-y-3">
                  <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider">
                    General
                  </h4>
                  <div className="space-y-2">
                    <ParamRow label="Control Mode" value={selectedController.control_mode} />
                    <ParamRow label="Rated Power" value={`${selectedController.rated_power_kw} kW`} />
                    <ParamRow label="Rated Speed" value={`${selectedController.rated_speed_rpm} RPM`} />
                  </div>
                </div>

                {/* Pitch */}
                <div className="space-y-3">
                  <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider">
                    Pitch Control
                  </h4>
                  <div className="space-y-2">
                    <ParamRow label="Min Pitch" value={`${selectedController.min_pitch_deg} deg`} />
                    <ParamRow label="Max Pitch" value={`${selectedController.max_pitch_deg} deg`} />
                    <ParamRow label="Max Pitch Rate" value={`${selectedController.max_pitch_rate_deg_s} deg/s`} />
                    <ParamRow label="Kp (Pitch)" value={selectedController.kp_pitch.toExponential(3)} />
                    <ParamRow label="Ki (Pitch)" value={selectedController.ki_pitch.toExponential(3)} />
                  </div>
                </div>

                {/* Torque */}
                <div className="space-y-3">
                  <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider">
                    Torque Control
                  </h4>
                  <div className="space-y-2">
                    <ParamRow label="Kp (Torque)" value={selectedController.kp_torque.toExponential(3)} />
                    <ParamRow label="Ki (Torque)" value={selectedController.ki_torque.toExponential(3)} />
                    <ParamRow label="Max Torque" value={`${selectedController.max_torque_nm.toLocaleString()} Nm`} />
                  </div>
                </div>

                {/* Generator */}
                <div className="space-y-3">
                  <h4 className="text-sm font-semibold text-slate-600 uppercase tracking-wider">
                    Generator
                  </h4>
                  <div className="space-y-2">
                    <ParamRow label="Min Gen Speed" value={`${selectedController.min_gen_speed_rpm} RPM`} />
                    <ParamRow label="Max Gen Speed" value={`${selectedController.max_gen_speed_rpm} RPM`} />
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

function ParamRow({ label, value }: { label: string; value: string }) {
  return (
    <div className="flex items-center justify-between rounded-md bg-slate-50 px-3 py-2">
      <span className="text-xs text-slate-500">{label}</span>
      <span className="text-sm font-medium text-slate-800 font-mono">
        {value}
      </span>
    </div>
  );
}
