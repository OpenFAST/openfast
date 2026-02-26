import { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { simulationsApi } from '@/api/client';
import type {
  Simulation,
  ResultsStatistics,
  ResultsDEL,
  ResultsExtreme,
} from '@/types';
import { BarChart3, Loader2, Download, Filter } from 'lucide-react';
import clsx from 'clsx';

interface ResultsData {
  statistics: ResultsStatistics[];
  dels: ResultsDEL[];
  extremes: ResultsExtreme[];
}

export default function ResultsDashboard() {
  const { projectId } = useParams<{ projectId: string }>();
  const [simulations, setSimulations] = useState<Simulation[]>([]);
  const [selectedSim, setSelectedSim] = useState<Simulation | null>(null);
  const [results, setResults] = useState<ResultsData | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [isLoadingResults, setIsLoadingResults] = useState(false);
  const [activeTab, setActiveTab] = useState<'statistics' | 'del' | 'extreme'>(
    'statistics',
  );

  useEffect(() => {
    if (!projectId) return;
    setIsLoading(true);
    simulationsApi
      .list(projectId)
      .then((data) => {
        const completed = data.filter((s) => s.status === 'completed');
        setSimulations(completed);
        if (completed.length > 0) setSelectedSim(completed[0]);
      })
      .catch(() => {})
      .finally(() => setIsLoading(false));
  }, [projectId]);

  useEffect(() => {
    if (!projectId || !selectedSim) return;
    setIsLoadingResults(true);
    simulationsApi
      .getResults(projectId, selectedSim.id)
      .then(setResults)
      .catch(() => setResults(null))
      .finally(() => setIsLoadingResults(false));
  }, [projectId, selectedSim]);

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
            Results Dashboard
          </h2>
          <p className="text-sm text-slate-500">
            Explore simulation statistics, DELs, and extreme loads
          </p>
        </div>
        {results && (
          <button className="btn-secondary">
            <Download className="h-4 w-4" />
            Export
          </button>
        )}
      </div>

      {simulations.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <BarChart3 className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">
            No completed simulations
          </h3>
          <p className="mt-1 text-sm text-slate-500 max-w-sm text-center">
            Run a simulation to completion, then view the results here.
          </p>
        </div>
      ) : (
        <div className="space-y-6">
          {/* Simulation selector */}
          <div className="flex items-center gap-4">
            <div className="flex items-center gap-2">
              <Filter className="h-4 w-4 text-slate-400" />
              <span className="text-sm font-medium text-slate-600">
                Simulation:
              </span>
            </div>
            <div className="flex gap-2">
              {simulations.map((sim) => (
                <button
                  key={sim.id}
                  onClick={() => setSelectedSim(sim)}
                  className={clsx(
                    'rounded-lg border px-4 py-2 text-sm font-medium transition-all',
                    selectedSim?.id === sim.id
                      ? 'border-accent-300 bg-accent-50 text-accent-700'
                      : 'border-slate-200 bg-white text-slate-600 hover:border-slate-300',
                  )}
                >
                  {sim.name}
                </button>
              ))}
            </div>
          </div>

          {isLoadingResults ? (
            <div className="flex items-center justify-center py-16">
              <Loader2 className="h-6 w-6 animate-spin text-accent-500" />
            </div>
          ) : results ? (
            <>
              {/* Tab toggle */}
              <div className="flex gap-1 rounded-lg bg-slate-100 p-1 max-w-md">
                {[
                  { key: 'statistics' as const, label: 'Statistics' },
                  { key: 'del' as const, label: 'Fatigue DELs' },
                  { key: 'extreme' as const, label: 'Extremes' },
                ].map((tab) => (
                  <button
                    key={tab.key}
                    onClick={() => setActiveTab(tab.key)}
                    className={clsx(
                      'flex-1 rounded-md px-3 py-1.5 text-sm font-medium transition-all',
                      activeTab === tab.key
                        ? 'bg-white text-slate-800 shadow-sm'
                        : 'text-slate-500 hover:text-slate-700',
                    )}
                  >
                    {tab.label}
                  </button>
                ))}
              </div>

              {/* Statistics */}
              {activeTab === 'statistics' &&
                results.statistics.length > 0 && (
                  <div className="rounded-xl border border-slate-200 bg-white overflow-hidden">
                    <div className="overflow-x-auto">
                      <table className="min-w-full text-sm">
                        <thead className="bg-slate-50">
                          <tr>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                              DLC
                            </th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                              Wind
                            </th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                              Channel
                            </th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">
                              Mean
                            </th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">
                              Std
                            </th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">
                              Min
                            </th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">
                              Max
                            </th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-slate-100">
                          {results.statistics.slice(0, 50).map((stat, si) =>
                            stat.statistics.map((ch, ci) => (
                              <tr
                                key={`${si}-${ci}`}
                                className="hover:bg-slate-50"
                              >
                                {ci === 0 && (
                                  <>
                                    <td
                                      className="px-3 py-2 font-medium text-slate-800"
                                      rowSpan={stat.statistics.length}
                                    >
                                      {stat.dlc_number}
                                    </td>
                                    <td
                                      className="px-3 py-2 text-slate-700"
                                      rowSpan={stat.statistics.length}
                                    >
                                      {stat.wind_speed_mps} m/s
                                    </td>
                                  </>
                                )}
                                <td className="px-3 py-2 text-slate-700 font-mono text-xs">
                                  {ch.channel}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {ch.mean.toExponential(3)}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {ch.std.toExponential(3)}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {ch.min.toExponential(3)}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {ch.max.toExponential(3)}
                                </td>
                              </tr>
                            )),
                          )}
                        </tbody>
                      </table>
                    </div>
                  </div>
                )}

              {/* DELs */}
              {activeTab === 'del' && results.dels.length > 0 && (
                <div className="rounded-xl border border-slate-200 bg-white overflow-hidden">
                  <div className="overflow-x-auto">
                    <table className="min-w-full text-sm">
                      <thead className="bg-slate-50">
                        <tr>
                          <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                            Channel
                          </th>
                          <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                            DLC
                          </th>
                          <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">
                            Wohler Exp.
                          </th>
                          <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">
                            Lifetime DEL
                          </th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-slate-100">
                        {results.dels.map((del, i) => (
                          <tr key={i} className="hover:bg-slate-50">
                            <td className="px-3 py-2 font-mono text-xs text-slate-800">
                              {del.channel}
                            </td>
                            <td className="px-3 py-2 text-slate-700">
                              {del.dlc_number}
                            </td>
                            <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                              {del.wohler_exponent}
                            </td>
                            <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs font-semibold">
                              {del.lifetime_del.toExponential(3)}
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
              )}

              {/* Extremes */}
              {activeTab === 'extreme' && results.extremes.length > 0 && (
                <div className="rounded-xl border border-slate-200 bg-white overflow-hidden">
                  <div className="overflow-x-auto">
                    <table className="min-w-full text-sm">
                      <thead className="bg-slate-50">
                        <tr>
                          <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                            Channel
                          </th>
                          <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                            DLC
                          </th>
                          <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">
                            Extreme Value
                          </th>
                          <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">
                            Time (s)
                          </th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-slate-100">
                        {results.extremes.map((ext, i) => (
                          <tr key={i} className="hover:bg-slate-50">
                            <td className="px-3 py-2 font-mono text-xs text-slate-800">
                              {ext.channel}
                            </td>
                            <td className="px-3 py-2 text-slate-700">
                              {ext.dlc_number}
                            </td>
                            <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs font-semibold">
                              {ext.extreme_value.toExponential(3)}
                            </td>
                            <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                              {ext.time_of_extreme.toFixed(2)}
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </div>
              )}

              {/* Empty tab states */}
              {activeTab === 'statistics' &&
                results.statistics.length === 0 && (
                  <EmptyTabState message="No statistics data available for this simulation." />
                )}
              {activeTab === 'del' && results.dels.length === 0 && (
                <EmptyTabState message="No fatigue DEL data available for this simulation." />
              )}
              {activeTab === 'extreme' &&
                results.extremes.length === 0 && (
                  <EmptyTabState message="No extreme load data available for this simulation." />
                )}
            </>
          ) : (
            <div className="rounded-xl border border-slate-200 bg-white p-8 text-center">
              <p className="text-sm text-slate-500">
                Select a completed simulation to view results.
              </p>
            </div>
          )}
        </div>
      )}
    </div>
  );
}

function EmptyTabState({ message }: { message: string }) {
  return (
    <div className="rounded-xl border border-slate-200 bg-white p-12 text-center">
      <BarChart3 className="mx-auto h-10 w-10 text-slate-300 mb-3" />
      <p className="text-sm text-slate-500">{message}</p>
    </div>
  );
}
