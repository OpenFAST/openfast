import { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { simulationsApi } from '@/api/client';
import type {
  Simulation,
  ResultsStatistics,
  ResultsDEL,
  ResultsExtreme,
} from '@/types';
import { BarChart3, Loader2 } from 'lucide-react';
import clsx from 'clsx';
import toast from 'react-hot-toast';

interface ResultsData {
  statistics: ResultsStatistics[];
  dels: ResultsDEL[];
  extremes: ResultsExtreme[];
}

export default function ResultsDashboard() {
  const { projectId } = useParams<{ projectId: string }>();
  const [simulations, setSimulations] = useState<Simulation[]>([]);
  const [selectedSimId, setSelectedSimId] = useState<string>('');
  const [results, setResults] = useState<ResultsData | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [isLoadingResults, setIsLoadingResults] = useState(false);
  const [activeTab, setActiveTab] = useState<'statistics' | 'del' | 'extreme'>('statistics');

  useEffect(() => {
    if (!projectId) return;
    setIsLoading(true);
    simulationsApi
      .list(projectId)
      .then((data) => {
        const completed = data.filter((s) => s.status === 'completed');
        setSimulations(completed);
        if (completed.length > 0) {
          setSelectedSimId(completed[0].id);
        }
      })
      .catch(() => toast.error('Failed to load simulations'))
      .finally(() => setIsLoading(false));
  }, [projectId]);

  useEffect(() => {
    if (!projectId || !selectedSimId) {
      setResults(null);
      return;
    }
    setIsLoadingResults(true);
    simulationsApi
      .getResults(projectId, selectedSimId)
      .then(setResults)
      .catch(() => {
        setResults(null);
        toast.error('Failed to load results');
      })
      .finally(() => setIsLoadingResults(false));
  }, [projectId, selectedSimId]);

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
          <h2 className="text-xl font-bold text-slate-800">Results Dashboard</h2>
          <p className="text-sm text-slate-500">
            Explore simulation statistics, DELs, and extreme loads
          </p>
        </div>
      </div>

      {simulations.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <BarChart3 className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">No completed simulations</h3>
          <p className="mt-1 text-sm text-slate-500 max-w-sm text-center">
            Run a simulation to completion, then view the results here.
          </p>
        </div>
      ) : (
        <div className="space-y-6">
          {/* Simulation selector */}
          <div className="flex items-center gap-4">
            <label className="text-sm font-medium text-slate-600">Simulation:</label>
            <select
              value={selectedSimId}
              onChange={(e) => setSelectedSimId(e.target.value)}
              className="rounded-lg border border-slate-300 px-4 py-2 text-sm focus:border-accent-500 focus:ring-1 focus:ring-accent-500"
            >
              {simulations.map((sim) => (
                <option key={sim.id} value={sim.id}>{sim.name}</option>
              ))}
            </select>
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
                  { key: 'extreme' as const, label: 'Extreme Loads' },
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

              {/* Statistics tab */}
              {activeTab === 'statistics' && (
                results.statistics.length > 0 ? (
                  <div className="rounded-xl border border-slate-200 bg-white overflow-hidden">
                    <div className="overflow-x-auto">
                      <table className="min-w-full text-sm">
                        <thead className="bg-slate-50">
                          <tr>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">DLC</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Wind (m/s)</th>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Channel</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">Mean</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">Std</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">Min</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">Max</th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-slate-100">
                          {results.statistics.map((stat) => {
                            if (!stat.channel_statistics) return null;
                            const entries = Object.entries(stat.channel_statistics);
                            return entries.map(([channel, vals], ci) => (
                              <tr key={`${stat.id}-${ci}`} className="hover:bg-slate-50">
                                {ci === 0 && (
                                  <>
                                    <td className="px-3 py-2 font-medium text-slate-800" rowSpan={entries.length}>
                                      {stat.dlc_number}
                                    </td>
                                    <td className="px-3 py-2 text-slate-700" rowSpan={entries.length}>
                                      {stat.wind_speed}
                                    </td>
                                  </>
                                )}
                                <td className="px-3 py-2 text-slate-700 font-mono text-xs">{channel}</td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {vals.mean.toExponential(3)}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {vals.std.toExponential(3)}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {vals.min.toExponential(3)}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {vals.max.toExponential(3)}
                                </td>
                              </tr>
                            ));
                          })}
                        </tbody>
                      </table>
                    </div>
                  </div>
                ) : (
                  <EmptyTabState message="No statistics data available for this simulation." />
                )
              )}

              {/* DEL tab */}
              {activeTab === 'del' && (
                results.dels.length > 0 ? (
                  <div className="rounded-xl border border-slate-200 bg-white overflow-hidden">
                    <div className="overflow-x-auto">
                      <table className="min-w-full text-sm">
                        <thead className="bg-slate-50">
                          <tr>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Channel</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">DEL Value</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">m Exponent</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">N Equivalent</th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-slate-100">
                          {results.dels.map((del) => {
                            if (!del.del_results) return null;
                            return Object.entries(del.del_results).map(([channel, value], i) => (
                              <tr key={`${del.id}-${i}`} className="hover:bg-slate-50">
                                <td className="px-3 py-2 font-mono text-xs text-slate-800">{channel}</td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs font-semibold">
                                  {value.toExponential(3)}
                                </td>
                                {i === 0 && (
                                  <>
                                    <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs" rowSpan={Object.keys(del.del_results!).length}>
                                      {del.m_exponent}
                                    </td>
                                    <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs" rowSpan={Object.keys(del.del_results!).length}>
                                      {del.n_equivalent.toExponential(2)}
                                    </td>
                                  </>
                                )}
                              </tr>
                            ));
                          })}
                        </tbody>
                      </table>
                    </div>
                  </div>
                ) : (
                  <EmptyTabState message="No fatigue DEL data available for this simulation." />
                )
              )}

              {/* Extreme tab */}
              {activeTab === 'extreme' && (
                results.extremes.length > 0 ? (
                  <div className="rounded-xl border border-slate-200 bg-white overflow-hidden">
                    <div className="overflow-x-auto">
                      <table className="min-w-full text-sm">
                        <thead className="bg-slate-50">
                          <tr>
                            <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">Channel</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">Max</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">Min</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">Safety Factor</th>
                            <th className="px-3 py-2 text-right text-xs font-medium text-slate-500 uppercase">Design Value</th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-slate-100">
                          {results.extremes.map((ext) => {
                            if (!ext.extreme_loads) return null;
                            return Object.entries(ext.extreme_loads).map(([channel, vals], i) => (
                              <tr key={`${ext.id}-${i}`} className="hover:bg-slate-50">
                                <td className="px-3 py-2 font-mono text-xs text-slate-800">{channel}</td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {vals.max.toExponential(3)}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {vals.min.toExponential(3)}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs">
                                  {vals.safety_factor.toFixed(2)}
                                </td>
                                <td className="px-3 py-2 text-right text-slate-700 font-mono text-xs font-semibold">
                                  {vals.design_value.toExponential(3)}
                                </td>
                              </tr>
                            ));
                          })}
                        </tbody>
                      </table>
                    </div>
                  </div>
                ) : (
                  <EmptyTabState message="No extreme load data available for this simulation." />
                )
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
