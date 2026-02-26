import { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { dlcDefinitionsApi } from '@/api/client';
import type { DLCDefinition } from '@/types';
import { Table2, Plus, Loader2, Trash2 } from 'lucide-react';
import clsx from 'clsx';

export default function DLCMatrix() {
  const { projectId } = useParams<{ projectId: string }>();
  const [dlcDefs, setDlcDefs] = useState<DLCDefinition[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedDlc, setSelectedDlc] = useState<DLCDefinition | null>(null);

  useEffect(() => {
    if (!projectId) return;
    setIsLoading(true);
    dlcDefinitionsApi
      .list(projectId)
      .then((data) => {
        setDlcDefs(data);
        if (data.length > 0) setSelectedDlc(data[0]);
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
          <h2 className="text-xl font-bold text-slate-800">DLC Matrix</h2>
          <p className="text-sm text-slate-500">
            Define Design Load Case sets for IEC certification
          </p>
        </div>
        <button className="btn-primary">
          <Plus className="h-4 w-4" />
          New DLC Set
        </button>
      </div>

      {dlcDefs.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Table2 className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">
            No DLC definitions
          </h3>
          <p className="mt-1 text-sm text-slate-500 max-w-sm text-center">
            Create Design Load Case definitions to configure the simulation
            matrix for IEC certification.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* DLC list */}
          <div className="space-y-2">
            {dlcDefs.map((dlc) => (
              <button
                key={dlc.id}
                onClick={() => setSelectedDlc(dlc)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedDlc?.id === dlc.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <span className="font-medium text-slate-800">{dlc.name}</span>
                <div className="mt-1 text-xs text-slate-500">
                  {dlc.iec_standard} &middot; {dlc.cases.length} case
                  {dlc.cases.length !== 1 ? 's' : ''}
                </div>
              </button>
            ))}
          </div>

          {/* DLC detail */}
          {selectedDlc && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-4">
                <div>
                  <h3 className="text-lg font-semibold text-slate-800">
                    {selectedDlc.name}
                  </h3>
                  <p className="text-sm text-slate-500">
                    Standard: {selectedDlc.iec_standard}
                  </p>
                </div>
                <button className="text-slate-400 hover:text-danger-500 transition-colors">
                  <Trash2 className="h-4 w-4" />
                </button>
              </div>

              {/* Cases table */}
              <div className="overflow-x-auto rounded-lg border border-slate-200">
                <table className="min-w-full text-sm">
                  <thead className="bg-slate-50">
                    <tr>
                      <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                        DLC
                      </th>
                      <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                        Wind (m/s)
                      </th>
                      <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                        Yaw Error
                      </th>
                      <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                        Seeds
                      </th>
                      <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                        Turb. Model
                      </th>
                      <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                        Fault
                      </th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-100">
                    {selectedDlc.cases.map((c, i) => (
                      <tr key={i} className="hover:bg-slate-50">
                        <td className="px-3 py-2 font-medium text-slate-800">
                          {c.dlc_number}
                        </td>
                        <td className="px-3 py-2 text-slate-700">
                          {c.wind_speed_mps}
                        </td>
                        <td className="px-3 py-2 text-slate-700">
                          {c.yaw_error_deg}&deg;
                        </td>
                        <td className="px-3 py-2 text-slate-700">
                          {c.num_seeds}
                        </td>
                        <td className="px-3 py-2 text-slate-700">
                          {c.turbsim_params.iec_turb_model}
                        </td>
                        <td className="px-3 py-2 text-slate-700">
                          {c.fault_type || '\u2014'}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>

              <div className="mt-4 rounded-lg bg-slate-50 p-3">
                <p className="text-xs text-slate-500">
                  Total simulations:{' '}
                  <span className="font-semibold text-slate-700">
                    {selectedDlc.cases.reduce(
                      (sum, c) => sum + c.num_seeds,
                      0,
                    )}
                  </span>{' '}
                  (across {selectedDlc.cases.length} case specifications)
                </p>
              </div>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
