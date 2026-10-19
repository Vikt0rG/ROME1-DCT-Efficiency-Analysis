#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include "utils.hpp"
#include "analysis/perFileHelpers.hpp"
#include "analysis/summaryHelpers.hpp"

namespace summaryHelpers {

void setupBranches(TTree* summary_tree, MeasurementMetadata& metadata, MeasurementData& data) {
    summary_tree->Branch("name", &metadata.name);
    summary_tree->Branch("measurement_type", &metadata.measurement_type);
    summary_tree->Branch("group_name", &metadata.group_name);
    summary_tree->Branch("mixture", &metadata.mixture);
    summary_tree->Branch("source", &metadata.source);
    summary_tree->Branch("filter", &metadata.filter);
    summary_tree->Branch("lv_setting", &metadata.lv_setting);
    summary_tree->Branch("scanned_layer", &metadata.scanned_layer);
    summary_tree->Branch("scanned_hv", &metadata.scanned_hv);
    summary_tree->Branch("other_hv", &metadata.other_hv);

    // Efficiency results
    summary_tree->Branch("eff_eta1_external_before_reco", &data.efficiency_results.eta1_efficiency_external, "eff_eta1_external_before_reco[3]/D");
    summary_tree->Branch("eff_eta2_external_before_reco", &data.efficiency_results.eta2_efficiency_external, "eff_eta2_external_before_reco[3]/D");
    summary_tree->Branch("eff_or_external_before_reco", &data.efficiency_results.eta_or_efficiency_external, "eff_or_external_before_reco[3]/D");
    summary_tree->Branch("eff_and_external_before_reco", &data.efficiency_results.eta_and_efficiency_external, "eff_and_external_before_reco[3]/D");
    summary_tree->Branch("eff_eta1_external_before_reco_error", &data.efficiency_results.eta1_efficiency_external_error, "eff_eta1_external_before_reco_error[6]/D");
    summary_tree->Branch("eff_eta2_external_before_reco_error", &data.efficiency_results.eta2_efficiency_external_error, "eff_eta2_external_before_reco_error[6]/D");
    summary_tree->Branch("eff_or_external_before_reco_error", &data.efficiency_results.eta_or_efficiency_external_error, "eff_or_external_before_reco_error[6]/D");
    summary_tree->Branch("eff_and_external_before_reco_error", &data.efficiency_results.eta_and_efficiency_external_error, "eff_and_external_before_reco_error[6]/D");

    summary_tree->Branch("eff_eta1_rpc_before_reco", &data.efficiency_results.eta1_efficiency_rpc, "eff_eta1_rpc_before_reco[3]/D");
    summary_tree->Branch("eff_eta2_rpc_before_reco", &data.efficiency_results.eta2_efficiency_rpc, "eff_eta2_rpc_before_reco[3]/D");
    summary_tree->Branch("eff_or_rpc_before_reco", &data.efficiency_results.eta_or_efficiency_rpc, "eff_or_rpc_before_reco[3]/D");
    summary_tree->Branch("eff_and_rpc_before_reco", &data.efficiency_results.eta_and_efficiency_rpc, "eff_and_rpc_before_reco[3]/D");
    summary_tree->Branch("eff_eta1_rpc_before_reco_error", &data.efficiency_results.eta1_efficiency_rpc_error, "eff_eta1_rpc_before_reco_error[6]/D");
    summary_tree->Branch("eff_eta2_rpc_before_reco_error", &data.efficiency_results.eta2_efficiency_rpc_error, "eff_eta2_rpc_before_reco_error[6]/D");
    summary_tree->Branch("eff_or_rpc_before_reco_error", &data.efficiency_results.eta_or_efficiency_rpc_error, "eff_or_rpc_before_reco_error[6]/D");
    summary_tree->Branch("eff_and_rpc_before_reco_error", &data.efficiency_results.eta_and_efficiency_rpc_error, "eff_and_rpc_before_reco_error[6]/D");

    summary_tree->Branch("eff_eta1_external_after_reco", &data.efficiency_results_tracks.eta1_efficiency_external, "eff_eta1_external_after_reco[3]/D");
    summary_tree->Branch("eff_eta2_external_after_reco", &data.efficiency_results_tracks.eta2_efficiency_external, "eff_eta2_external_after_reco[3]/D");
    summary_tree->Branch("eff_or_external_after_reco", &data.efficiency_results_tracks.eta_or_efficiency_external, "eff_or_external_after_reco[3]/D");
    summary_tree->Branch("eff_and_external_after_reco", &data.efficiency_results_tracks.eta_and_efficiency_external, "eff_and_external_after_reco[3]/D");
    summary_tree->Branch("eff_eta1_external_after_reco_error", &data.efficiency_results_tracks.eta1_efficiency_external_error, "eff_eta1_external_after_reco_error[6]/D");
    summary_tree->Branch("eff_eta2_external_after_reco_error", &data.efficiency_results_tracks.eta2_efficiency_external_error, "eff_eta2_external_after_reco_error[6]/D");
    summary_tree->Branch("eff_or_external_after_reco_error", &data.efficiency_results_tracks.eta_or_efficiency_external_error, "eff_or_external_after_reco_error[6]/D");
    summary_tree->Branch("eff_and_external_after_reco_error", &data.efficiency_results_tracks.eta_and_efficiency_external_error, "eff_and_external_after_reco_error[6]/D");

    summary_tree->Branch("eff_eta1_rpc_after_reco", &data.efficiency_results_tracks.eta1_efficiency_rpc, "eff_eta1_rpc_after_reco[3]/D");
    summary_tree->Branch("eff_eta2_rpc_after_reco", &data.efficiency_results_tracks.eta2_efficiency_rpc, "eff_eta2_rpc_after_reco[3]/D");
    summary_tree->Branch("eff_or_rpc_after_reco", &data.efficiency_results_tracks.eta_or_efficiency_rpc, "eff_or_rpc_after_reco[3]/D");
    summary_tree->Branch("eff_and_rpc_after_reco", &data.efficiency_results_tracks.eta_and_efficiency_rpc, "eff_and_rpc_after_reco[3]/D");
    summary_tree->Branch("eff_eta1_rpc_after_reco_error", &data.efficiency_results_tracks.eta1_efficiency_rpc_error, "eff_eta1_rpc_after_reco_error[6]/D");
    summary_tree->Branch("eff_eta2_rpc_after_reco_error", &data.efficiency_results_tracks.eta2_efficiency_rpc_error, "eff_eta2_rpc_after_reco_error[6]/D");
    summary_tree->Branch("eff_or_rpc_after_reco_error", &data.efficiency_results_tracks.eta_or_efficiency_rpc_error, "eff_or_rpc_after_reco_error[6]/D");
    summary_tree->Branch("eff_and_rpc_after_reco_error", &data.efficiency_results_tracks.eta_and_efficiency_rpc_error, "eff_and_rpc_after_reco_error[6]/D");

    // Cluster size results
    summary_tree->Branch("avg_cluster_size_eta1", &data.cluster_size_results.avg_cluster_size_eta1);
    summary_tree->Branch("avg_cluster_size_eta2", &data.cluster_size_results.avg_cluster_size_eta2);
    summary_tree->Branch("avg_cluster_size_eta1_layers", &data.cluster_size_results.avg_cluster_size_eta1_layers, "avg_cluster_size_eta1_layers[3]/D");
    summary_tree->Branch("avg_cluster_size_eta2_layers", &data.cluster_size_results.avg_cluster_size_eta2_layers, "avg_cluster_size_eta2_layers[3]/D");
    summary_tree->Branch("avg_cluster_size_eta1_error", &data.cluster_size_results.avg_cluster_size_eta1_error, "avg_cluster_size_eta1_error[2]/D");
    summary_tree->Branch("avg_cluster_size_eta2_error", &data.cluster_size_results.avg_cluster_size_eta2_error, "avg_cluster_size_eta2_error[2]/D");
    summary_tree->Branch("avg_cluster_size_eta1_layers_error", &data.cluster_size_results.avg_cluster_size_eta1_layers_error, "avg_cluster_size_eta1_layers_error[6]/D");
    summary_tree->Branch("avg_cluster_size_eta2_layers_error", &data.cluster_size_results.avg_cluster_size_eta2_layers_error, "avg_cluster_size_eta2_layers_error[6]/D");

    // Rate results
    summary_tree->Branch("rate_eta1_before_reco", &data.rate_results.rate_eta1);
    summary_tree->Branch("rate_eta2_before_reco", &data.rate_results.rate_eta2);
    summary_tree->Branch("rate_or_before_reco", &data.rate_results.rate_or);
    summary_tree->Branch("rate_and_before_reco", &data.rate_results.rate_and);
    summary_tree->Branch("rate_eta1_before_reco_error", &data.rate_results.rate_eta1_error, Form("rate_eta1_before_reco_error[%d]/D", 2));
    summary_tree->Branch("rate_eta2_before_reco_error", &data.rate_results.rate_eta2_error, Form("rate_eta2_before_reco_error[%d]/D", 2));
    summary_tree->Branch("rate_or_before_reco_error", &data.rate_results.rate_or_error, Form("rate_or_before_reco_error[%d]/D", 2));
    summary_tree->Branch("rate_and_before_reco_error", &data.rate_results.rate_and_error, Form("rate_and_before_reco_error[%d]/D", 2));

    summary_tree->Branch("rate_eta1_after_reco", &data.rate_results_tracks.rate_eta1);
    summary_tree->Branch("rate_eta2_after_reco", &data.rate_results_tracks.rate_eta2);
    summary_tree->Branch("rate_or_after_reco", &data.rate_results_tracks.rate_or);
    summary_tree->Branch("rate_and_after_reco", &data.rate_results_tracks.rate_and);
    summary_tree->Branch("rate_eta1_after_reco_error", &data.rate_results_tracks.rate_eta1_error, Form("rate_eta1_after_reco_error[%d]/D", 2));
    summary_tree->Branch("rate_eta2_after_reco_error", &data.rate_results_tracks.rate_eta2_error, Form("rate_eta2_after_reco_error[%d]/D", 2));
    summary_tree->Branch("rate_or_after_reco_error", &data.rate_results_tracks.rate_or_error, Form("rate_or_after_reco_error[%d]/D", 2));
    summary_tree->Branch("rate_and_after_reco_error", &data.rate_results_tracks.rate_and_error, Form("rate_and_after_reco_error[%d]/D", 2));

    summary_tree->Branch("rate_layers_eta1_before_reco", &data.rate_results.rate_layers_eta1, Form("rate_layers_eta1_before_reco[%d]/D", 3));
    summary_tree->Branch("rate_layers_eta2_before_reco", &data.rate_results.rate_layers_eta2, Form("rate_layers_eta2_before_reco[%d]/D", 3));
    summary_tree->Branch("rate_layers_or_before_reco", &data.rate_results.rate_layers_or, Form("rate_layers_or_before_reco[%d]/D", 3));
    summary_tree->Branch("rate_layers_and_before_reco", &data.rate_results.rate_layers_and, Form("rate_layers_and_before_reco[%d]/D", 3));
    summary_tree->Branch("rate_layers_eta1_before_reco_error", &data.rate_results.rate_layers_eta1_error, Form("rate_layers_eta1_before_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("rate_layers_eta2_before_reco_error", &data.rate_results.rate_layers_eta2_error, Form("rate_layers_eta2_before_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("rate_layers_or_before_reco_error", &data.rate_results.rate_layers_or_error, Form("rate_layers_or_before_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("rate_layers_and_before_reco_error", &data.rate_results.rate_layers_and_error, Form("rate_layers_and_before_reco_error[%d]/D", LAYER_COUNT * 2));

    summary_tree->Branch("rate_layers_eta1_after_reco", &data.rate_results_tracks.rate_layers_eta1, Form("rate_layers_eta1_after_reco[%d]/D", 3));
    summary_tree->Branch("rate_layers_eta2_after_reco", &data.rate_results_tracks.rate_layers_eta2, Form("rate_layers_eta2_after_reco[%d]/D", 3));
    summary_tree->Branch("rate_layers_or_after_reco", &data.rate_results_tracks.rate_layers_or, Form("rate_layers_or_after_reco[%d]/D", 3));
    summary_tree->Branch("rate_layers_and_after_reco", &data.rate_results_tracks.rate_layers_and, Form("rate_layers_and_after_reco[%d]/D", 3));
    summary_tree->Branch("rate_layers_eta1_after_reco_error", &data.rate_results_tracks.rate_layers_eta1_error, Form("rate_layers_eta1_after_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("rate_layers_eta2_after_reco_error", &data.rate_results_tracks.rate_layers_eta2_error, Form("rate_layers_eta2_after_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("rate_layers_or_after_reco_error", &data.rate_results_tracks.rate_layers_or_error, Form("rate_layers_or_after_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("rate_layers_and_after_reco_error", &data.rate_results_tracks.rate_layers_and_error, Form("rate_layers_and_after_reco_error[%d]/D", LAYER_COUNT * 2));

    summary_tree->Branch("rate_strips_eta1_before_reco", &data.rate_results.rate_strips_eta1, Form("rate_strips_eta1_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("rate_strips_eta2_before_reco", &data.rate_results.rate_strips_eta2, Form("rate_strips_eta2_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("rate_strips_or_before_reco", &data.rate_results.rate_strips_or, Form("rate_strips_or_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("rate_strips_and_before_reco", &data.rate_results.rate_strips_and, Form("rate_strips_and_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("rate_strips_eta1_before_reco_error", &data.rate_results.rate_strips_eta1_error, Form("rate_strips_eta1_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("rate_strips_eta2_before_reco_error", &data.rate_results.rate_strips_eta2_error, Form("rate_strips_eta2_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("rate_strips_or_before_reco_error", &data.rate_results.rate_strips_or_error, Form("rate_strips_or_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("rate_strips_and_before_reco_error", &data.rate_results.rate_strips_and_error, Form("rate_strips_and_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("rate_strips_eta1_after_reco", &data.rate_results_tracks.rate_strips_eta1, Form("rate_strips_eta1_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("rate_strips_eta2_after_reco", &data.rate_results_tracks.rate_strips_eta2, Form("rate_strips_eta2_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("rate_strips_or_after_reco", &data.rate_results_tracks.rate_strips_or, Form("rate_strips_or_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("rate_strips_and_after_reco", &data.rate_results_tracks.rate_strips_and, Form("rate_strips_and_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("rate_strips_eta1_after_reco_error", &data.rate_results_tracks.rate_strips_eta1_error, Form("rate_strips_eta1_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("rate_strips_eta2_after_reco_error", &data.rate_results_tracks.rate_strips_eta2_error, Form("rate_strips_eta2_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("rate_strips_or_after_reco_error", &data.rate_results_tracks.rate_strips_or_error, Form("rate_strips_or_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("rate_strips_and_after_reco_error", &data.rate_results_tracks.rate_strips_and_error, Form("rate_strips_and_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    // ToT results
    summary_tree->Branch("avg_tot_layer_eta1_before_reco", &data.tot_results.avg_tot_eta1_layer, Form("avg_tot_layer_eta1_before_reco[%d]/D", LAYER_COUNT));
    summary_tree->Branch("avg_tot_layer_eta2_before_reco", &data.tot_results.avg_tot_eta2_layer, Form("avg_tot_layer_eta2_before_reco[%d]/D", LAYER_COUNT));
    summary_tree->Branch("avg_tot_layer_eta1_before_reco_error", &data.tot_results.avg_tot_eta1_layer_error, Form("avg_tot_layer_eta1_before_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("avg_tot_layer_eta2_before_reco_error", &data.tot_results.avg_tot_eta2_layer_error, Form("avg_tot_layer_eta2_before_reco_error[%d]/D", LAYER_COUNT * 2));

    summary_tree->Branch("avg_tot_layer_eta1_after_reco", &data.tot_results_tracks.avg_tot_eta1_layer, Form("avg_tot_layer_eta1_after_reco[%d]/D", LAYER_COUNT));
    summary_tree->Branch("avg_tot_layer_eta2_after_reco", &data.tot_results_tracks.avg_tot_eta2_layer, Form("avg_tot_layer_eta2_after_reco[%d]/D", LAYER_COUNT));
    summary_tree->Branch("avg_tot_layer_eta1_after_reco_error", &data.tot_results_tracks.avg_tot_eta1_layer_error, Form("avg_tot_layer_eta1_after_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("avg_tot_layer_eta2_after_reco_error", &data.tot_results_tracks.avg_tot_eta2_layer_error, Form("avg_tot_layer_eta2_after_reco_error[%d]/D", LAYER_COUNT * 2));

    summary_tree->Branch("beam_region_avg_tot_layer_eta1_before_reco", &data.tot_results_beam.avg_tot_eta1_layer, Form("beam_region_avg_tot_layer_eta1_before_reco[%d]/D", LAYER_COUNT));
    summary_tree->Branch("beam_region_avg_tot_layer_eta2_before_reco", &data.tot_results_beam.avg_tot_eta2_layer, Form("beam_region_avg_tot_layer_eta2_before_reco[%d]/D", LAYER_COUNT));
    summary_tree->Branch("beam_region_avg_tot_layer_eta1_before_reco_error", &data.tot_results_beam.avg_tot_eta1_layer_error, Form("beam_region_avg_tot_layer_eta1_before_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("beam_region_avg_tot_layer_eta2_before_reco_error", &data.tot_results_beam.avg_tot_eta2_layer_error, Form("beam_region_avg_tot_layer_eta2_before_reco_error[%d]/D", LAYER_COUNT * 2));

    summary_tree->Branch("beam_region_avg_tot_layer_eta1_after_reco", &data.tot_results_beam_track.avg_tot_eta1_layer, Form("beam_region_avg_tot_layer_eta1_after_reco[%d]/D", LAYER_COUNT));
    summary_tree->Branch("beam_region_avg_tot_layer_eta2_after_reco", &data.tot_results_beam_track.avg_tot_eta2_layer, Form("beam_region_avg_tot_layer_eta2_after_reco[%d]/D", LAYER_COUNT));
    summary_tree->Branch("beam_region_avg_tot_layer_eta1_after_reco_error", &data.tot_results_beam_track.avg_tot_eta1_layer_error, Form("beam_region_avg_tot_layer_eta1_after_reco_error[%d]/D", LAYER_COUNT * 2));
    summary_tree->Branch("beam_region_avg_tot_layer_eta2_after_reco_error", &data.tot_results_beam_track.avg_tot_eta2_layer_error, Form("beam_region_avg_tot_layer_eta2_after_reco_error[%d]/D", LAYER_COUNT * 2));

    summary_tree->Branch("avg_tot_strip_eta1_before_reco", &data.tot_results.avg_tot_eta1, Form("avg_tot_strip_eta1_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_tot_strip_eta2_before_reco", &data.tot_results.avg_tot_eta2, Form("avg_tot_strip_eta2_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_tot_strip_eta1_before_reco_error", &data.tot_results.avg_tot_eta1_error, Form("avg_tot_strip_eta1_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("avg_tot_strip_eta2_before_reco_error", &data.tot_results.avg_tot_eta2_error, Form("avg_tot_strip_eta2_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("avg_tot_strip_eta1_after_reco", &data.tot_results_tracks.avg_tot_eta1, Form("avg_tot_strip_eta1_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_tot_strip_eta2_after_reco", &data.tot_results_tracks.avg_tot_eta2, Form("avg_tot_strip_eta2_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_tot_strip_eta1_after_reco_error", &data.tot_results_tracks.avg_tot_eta1_error, Form("avg_tot_strip_eta1_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("avg_tot_strip_eta2_after_reco_error", &data.tot_results_tracks.avg_tot_eta2_error, Form("avg_tot_strip_eta2_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("beam_region_avg_tot_strip_eta1_before_reco", &data.tot_results_beam.avg_tot_eta1, Form("beam_region_avg_tot_strip_eta1_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("beam_region_avg_tot_strip_eta2_before_reco", &data.tot_results_beam.avg_tot_eta2, Form("beam_region_avg_tot_strip_eta2_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("beam_region_avg_tot_strip_eta1_before_reco_error", &data.tot_results_beam.avg_tot_eta1_error, Form("beam_region_avg_tot_strip_eta1_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("beam_region_avg_tot_strip_eta2_before_reco_error", &data.tot_results_beam.avg_tot_eta2_error, Form("beam_region_avg_tot_strip_eta2_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("beam_region_avg_tot_strip_eta1_after_reco", &data.tot_results_beam_track.avg_tot_eta1, Form("beam_region_avg_tot_strip_eta1_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("beam_region_avg_tot_strip_eta2_after_reco", &data.tot_results_beam_track.avg_tot_eta2, Form("beam_region_avg_tot_strip_eta2_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("beam_region_avg_tot_strip_eta1_after_reco_error", &data.tot_results_beam_track.avg_tot_eta1_error, Form("beam_region_avg_tot_strip_eta1_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("beam_region_avg_tot_strip_eta2_after_reco_error", &data.tot_results_beam_track.avg_tot_eta2_error, Form("beam_region_avg_tot_strip_eta2_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    // Multiplicity results
    summary_tree->Branch("avg_multiplicity_eta1_before_reco", &data.multiplicity_results.avg_multiplicity_eta1, Form("avg_multiplicity_eta1_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_multiplicity_eta2_before_reco", &data.multiplicity_results.avg_multiplicity_eta2, Form("avg_multiplicity_eta2_before_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_multiplicity_eta1_before_reco_error", &data.multiplicity_results.avg_multiplicity_eta1_error, Form("avg_multiplicity_eta1_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("avg_multiplicity_eta2_before_reco_error", &data.multiplicity_results.avg_multiplicity_eta2_error, Form("avg_multiplicity_eta2_before_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    summary_tree->Branch("avg_multiplicity_eta1_after_reco", &data.multiplicity_results_tracks.avg_multiplicity_eta1, Form("avg_multiplicity_eta1_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_multiplicity_eta2_after_reco", &data.multiplicity_results_tracks.avg_multiplicity_eta2, Form("avg_multiplicity_eta2_after_reco[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER));
    summary_tree->Branch("avg_multiplicity_eta1_after_reco_error", &data.multiplicity_results_tracks.avg_multiplicity_eta1_error, Form("avg_multiplicity_eta1_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));
    summary_tree->Branch("avg_multiplicity_eta2_after_reco_error", &data.multiplicity_results_tracks.avg_multiplicity_eta2_error, Form("avg_multiplicity_eta2_after_reco_error[%d][%d]/D", LAYER_COUNT, STRIPS_PER_LAYER * 2));

    // ToF results
    for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {
        summary_tree->Branch(Form("time_of_flight_layer_%s_eta1", LAYER_PAIR_SUFFIXES[i].c_str()), &data.tof_results.time_of_flight_eta1[i]);
        summary_tree->Branch(Form("time_of_flight_layer_%s_eta2", LAYER_PAIR_SUFFIXES[i].c_str()), &data.tof_results.time_of_flight_eta2[i]);

        summary_tree->Branch(Form("avg_time_of_flight_layer_%s_eta1", LAYER_PAIR_SUFFIXES[i].c_str()), &data.tof_results.avg_time_of_flight_eta1[i]);
        summary_tree->Branch(Form("avg_time_of_flight_layer_%s_eta2", LAYER_PAIR_SUFFIXES[i].c_str()), &data.tof_results.avg_time_of_flight_eta2[i]);

        summary_tree->Branch(Form("avg_time_of_flight_layer_%s_eta1_error", LAYER_PAIR_SUFFIXES[i].c_str()), &data.tof_results.avg_time_of_flight_eta1_error[i], Form("avg_time_of_flight_layer_%s_eta1_error[%d]/D", LAYER_PAIR_SUFFIXES[i].c_str(), 2));
        summary_tree->Branch(Form("avg_time_of_flight_layer_%s_eta2_error", LAYER_PAIR_SUFFIXES[i].c_str()), &data.tof_results.avg_time_of_flight_eta2_error[i], Form("avg_time_of_flight_layer_%s_eta2_error[%d]/D", LAYER_PAIR_SUFFIXES[i].c_str(), 2));

        summary_tree->Branch(Form("time_resolution_layer_%s_eta1", LAYER_PAIR_SUFFIXES[i].c_str()), &data.time_resolution_results.time_resolution_global_eta1[i]);
        summary_tree->Branch(Form("time_resolution_layer_%s_eta2", LAYER_PAIR_SUFFIXES[i].c_str()), &data.time_resolution_results.time_resolution_global_eta2[i]);
        summary_tree->Branch(Form("time_resolution_layer_%s_eta1_error", LAYER_PAIR_SUFFIXES[i].c_str()), &data.time_resolution_results.time_resolution_global_eta1_error[i], Form("time_resolution_layer_%s_eta1_error[%d]/D", LAYER_PAIR_SUFFIXES[i].c_str(), 2));
        summary_tree->Branch(Form("time_resolution_layer_%s_eta2_error", LAYER_PAIR_SUFFIXES[i].c_str()), &data.time_resolution_results.time_resolution_global_eta2_error[i], Form("time_resolution_layer_%s_eta2_error[%d]/D", LAYER_PAIR_SUFFIXES[i].c_str(), 2));

        summary_tree->Branch(Form("time_resolution_layer_%s_strip_eta1", LAYER_PAIR_SUFFIXES[i].c_str()), &data.time_resolution_results.time_resolution_strip_eta1[i], Form("time_resolution_layer_%s_strip_eta1[%d]/D", LAYER_PAIR_SUFFIXES[i].c_str(), STRIPS_PER_LAYER));
        summary_tree->Branch(Form("time_resolution_layer_%s_strip_eta2", LAYER_PAIR_SUFFIXES[i].c_str()), &data.time_resolution_results.time_resolution_strip_eta2[i], Form("time_resolution_layer_%s_strip_eta2[%d]/D", LAYER_PAIR_SUFFIXES[i].c_str(), STRIPS_PER_LAYER));
        summary_tree->Branch(Form("time_resolution_layer_%s_strip_eta1_error", LAYER_PAIR_SUFFIXES[i].c_str()), &data.time_resolution_results.time_resolution_strip_eta1_error[i], Form("time_resolution_layer_%s_strip_eta1_error[%d]/D", LAYER_PAIR_SUFFIXES[i].c_str(), STRIPS_PER_LAYER * 2));
        summary_tree->Branch(Form("time_resolution_layer_%s_strip_eta2_error", LAYER_PAIR_SUFFIXES[i].c_str()), &data.time_resolution_results.time_resolution_strip_eta2_error[i], Form("time_resolution_layer_%s_strip_eta2_error[%d]/D", LAYER_PAIR_SUFFIXES[i].c_str(), STRIPS_PER_LAYER * 2));
    }
}

void getEfficiency(TFile* input_file, EfficiencyResults& eff_results, EfficiencyResults& track_eff_results) {

    auto extractGraphData = [](TGraphAsymmErrors* graph, double effs[4], ErrorRange errs[4]) {
        for (int i = 0; i < 4; ++i) {
            double dummy_x = 0.0;
            // GetPoint stores the Y value directly inside the target efficiency index array
            graph->GetPoint(i, dummy_x, effs[i]);

            // Extract the separate low and high bounds for the ErrorRange struct
            errs[i].low  = graph->GetErrorYlow(i);
            errs[i].high = graph->GetErrorYhigh(i);
        }
    };

    for (int layer = 0; layer < LAYER_COUNT; ++layer) {
        // ------------------------------------------------------------------------------
        // 1. External trigger only
        std::string path_ext = "efficiency_graphs/external_trigger/eff_external_trigger_layer" + std::to_string(layer);
        if (auto graph_ext = input_file->Get<TGraphAsymmErrors>(path_ext.c_str())) {
            double effs[4]; ErrorRange errs[4];
            extractGraphData(graph_ext, effs, errs);

            eff_results.eta1_efficiency_external[layer] = effs[0];
            eff_results.eta2_efficiency_external[layer] = effs[1];
            eff_results.eta_or_efficiency_external[layer] = effs[2];
            eff_results.eta_and_efficiency_external[layer] = effs[3];

            eff_results.eta1_efficiency_external_error[layer] = errs[0];
            eff_results.eta2_efficiency_external_error[layer] = errs[1];
            eff_results.eta_or_efficiency_external_error[layer] = errs[2];
            eff_results.eta_and_efficiency_external_error[layer] = errs[3];
        }

        // ------------------------------------------------------------------------------
        // 2. External trigger + RPC
        std::string path_rpc = "efficiency_graphs/external_plus_rpc_trigger/eff_rpc_layer" + std::to_string(layer);
        if (auto graph_rpc = input_file->Get<TGraphAsymmErrors>(path_rpc.c_str())) {
            double effs[4]; ErrorRange errs[4];
            extractGraphData(graph_rpc, effs, errs);

            eff_results.eta1_efficiency_rpc[layer] = effs[0];
            eff_results.eta2_efficiency_rpc[layer] = effs[1];
            eff_results.eta_or_efficiency_rpc[layer] = effs[2];
            eff_results.eta_and_efficiency_rpc[layer] = effs[3];

            eff_results.eta1_efficiency_rpc_error[layer] = errs[0];
            eff_results.eta2_efficiency_rpc_error[layer] = errs[1];
            eff_results.eta_or_efficiency_rpc_error[layer] = errs[2];
            eff_results.eta_and_efficiency_rpc_error[layer] = errs[3];
        }

        // ------------------------------------------------------------------------------
        // 3. Track external trigger only
        std::string path_track = "efficiency_graphs/track_external_trigger/track_eff_external_trigger_layer" + std::to_string(layer);
        if (auto graph_track = input_file->Get<TGraphAsymmErrors>(path_track.c_str())) {
            double effs[4]; ErrorRange errs[4];
            extractGraphData(graph_track, effs, errs);

            track_eff_results.eta1_efficiency_external[layer] = effs[0];
            track_eff_results.eta2_efficiency_external[layer] = effs[1];
            track_eff_results.eta_or_efficiency_external[layer] = effs[2];
            track_eff_results.eta_and_efficiency_external[layer] = effs[3];

            track_eff_results.eta1_efficiency_external_error[layer] = errs[0];
            track_eff_results.eta2_efficiency_external_error[layer] = errs[1];
            track_eff_results.eta_or_efficiency_external_error[layer] = errs[2];
            track_eff_results.eta_and_efficiency_external_error[layer] = errs[3];
        }

        // ------------------------------------------------------------------------------
        // 4. Track external trigger + RPC
        std::string path_track_rpc = "efficiency_graphs/track_external_plus_rpc_trigger/track_eff_rpc_layer" + std::to_string(layer);
        if (auto graph_track_rpc = input_file->Get<TGraphAsymmErrors>(path_track_rpc.c_str())) {
            double effs[4]; ErrorRange errs[4];
            extractGraphData(graph_track_rpc, effs, errs);

            track_eff_results.eta1_efficiency_rpc[layer] = effs[0];
            track_eff_results.eta2_efficiency_rpc[layer] = effs[1];
            track_eff_results.eta_or_efficiency_rpc[layer] = effs[2];
            track_eff_results.eta_and_efficiency_rpc[layer] = effs[3];

            track_eff_results.eta1_efficiency_rpc_error[layer] = errs[0];
            track_eff_results.eta2_efficiency_rpc_error[layer] = errs[1];
            track_eff_results.eta_or_efficiency_rpc_error[layer] = errs[2];
            track_eff_results.eta_and_efficiency_rpc_error[layer] = errs[3];
        }
    }
}

void getClusterSize(TFile* input_file, ClusterSizeResults& cluster_results) {
    TTree* cluster_tree = input_file->Get<TTree>("Clusterization");
    if (!cluster_tree) {
        std::cerr << "Error: Clusterization tree not found in file " << input_file->GetName() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    TTreeReader reader(cluster_tree);
    TTreeReaderValue<std::vector<int>> cluster_eta1(reader, "cluster_size_eta1");
    TTreeReaderValue<std::vector<int>> cluster_eta2(reader, "cluster_size_eta2");
    TTreeReaderValue<std::vector<int>> cluster_eta1_layer0(reader, "cluster_size_eta1_layer0");
    TTreeReaderValue<std::vector<int>> cluster_eta1_layer1(reader, "cluster_size_eta1_layer1");
    TTreeReaderValue<std::vector<int>> cluster_eta1_layer2(reader, "cluster_size_eta1_layer2");
    TTreeReaderValue<std::vector<int>> cluster_eta2_layer0(reader, "cluster_size_eta2_layer0");
    TTreeReaderValue<std::vector<int>> cluster_eta2_layer1(reader, "cluster_size_eta2_layer1");
    TTreeReaderValue<std::vector<int>> cluster_eta2_layer2(reader, "cluster_size_eta2_layer2");

    long long total_clusters_eta1 = 0, total_clusters_eta2 = 0;
    long long total_cluster_size_eta1 = 0, total_cluster_size_eta2 = 0;
    long long sum_sq_cluster_size_eta1 = 0, sum_sq_cluster_size_eta2 = 0;

    long long total_clusters_eta1_layers[3] = {0}, total_clusters_eta2_layers[3] = {0};
    long long total_cluster_size_eta1_layers[3] = {0}, total_cluster_size_eta2_layers[3] = {0};
    long long sum_sq_cluster_size_eta1_layers[3] = {0}, sum_sq_cluster_size_eta2_layers[3] = {0};

    auto accumulateBranch = [](auto& reader_val, long long& count, long long& sum, long long& sum_sq) {
        if (reader_val.GetSetupStatus() == 0) {
            count += static_cast<long long>(reader_val->size());
            for (int size : *reader_val) {
                sum += size;
                sum_sq += (static_cast<long long>(size) * size);
            }
        }
    };

    while (reader.Next()) {
        accumulateBranch(cluster_eta1, total_clusters_eta1, total_cluster_size_eta1, sum_sq_cluster_size_eta1);
        accumulateBranch(cluster_eta2, total_clusters_eta2, total_cluster_size_eta2, sum_sq_cluster_size_eta2);

        accumulateBranch(cluster_eta1_layer0, total_clusters_eta1_layers[0], total_cluster_size_eta1_layers[0], sum_sq_cluster_size_eta1_layers[0]);
        accumulateBranch(cluster_eta1_layer1, total_clusters_eta1_layers[1], total_cluster_size_eta1_layers[1], sum_sq_cluster_size_eta1_layers[1]);
        accumulateBranch(cluster_eta1_layer2, total_clusters_eta1_layers[2], total_cluster_size_eta1_layers[2], sum_sq_cluster_size_eta1_layers[2]);

        accumulateBranch(cluster_eta2_layer0, total_clusters_eta2_layers[0], total_cluster_size_eta2_layers[0], sum_sq_cluster_size_eta2_layers[0]);
        accumulateBranch(cluster_eta2_layer1, total_clusters_eta2_layers[1], total_cluster_size_eta2_layers[1], sum_sq_cluster_size_eta2_layers[1]);
        accumulateBranch(cluster_eta2_layer2, total_clusters_eta2_layers[2], total_cluster_size_eta2_layers[2], sum_sq_cluster_size_eta2_layers[2]);
    }

    auto calcMeanAndError = [](long long N, long long sum, long long sum_sq, double& mean_out, ErrorRange& err_out) {
        if (N <= 0) {
            mean_out = 0.0;
            err_out = ErrorRange{0.0};
            return;
        }
        mean_out = static_cast<double>(sum) / N;
        if (N > 1) {
            double variance = (static_cast<double>(sum_sq) - (static_cast<double>(sum) * sum) / N) / (N - 1);
            if (variance < 0.0) variance = 0.0;
            double std_dev = std::sqrt(variance);
            err_out = ErrorRange{std_dev / std::sqrt(N)};
        } else {
            err_out = ErrorRange{0.0};
        }
    };

    calcMeanAndError(total_clusters_eta1, total_cluster_size_eta1, sum_sq_cluster_size_eta1,
                        cluster_results.avg_cluster_size_eta1, cluster_results.avg_cluster_size_eta1_error);
    calcMeanAndError(total_clusters_eta2, total_cluster_size_eta2, sum_sq_cluster_size_eta2,
                        cluster_results.avg_cluster_size_eta2, cluster_results.avg_cluster_size_eta2_error);

    for (int layer_idx = 0; layer_idx < 3; ++layer_idx) {
        calcMeanAndError(total_clusters_eta1_layers[layer_idx], total_cluster_size_eta1_layers[layer_idx],
                            sum_sq_cluster_size_eta1_layers[layer_idx], cluster_results.avg_cluster_size_eta1_layers[layer_idx],
                            cluster_results.avg_cluster_size_eta1_layers_error[layer_idx]);
        calcMeanAndError(total_clusters_eta2_layers[layer_idx], total_cluster_size_eta2_layers[layer_idx],
                            sum_sq_cluster_size_eta2_layers[layer_idx], cluster_results.avg_cluster_size_eta2_layers[layer_idx],
                            cluster_results.avg_cluster_size_eta2_layers_error[layer_idx]);
    }
}

void getRate(TFile* input_file, RateResults& rate_results, bool in_valid_track_only) {
    TTree* processed_tree = dynamic_cast<TTree*>(input_file->Get("ProcessedData"));
    TTree* track_tree = dynamic_cast<TTree*>(input_file->Get("TrackReconstruction"));
    if (!processed_tree || !track_tree) {
        std::cerr << "Error: One or more required trees not found in file " << input_file->GetName() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    TTreeReader readerProcData(processed_tree);
    TTreeReader readerTrackData(track_tree);

    TTreeReaderValue<std::vector<int>> strips(readerProcData, "proc_strip");
    TTreeReaderValue<std::vector<int>> layers(readerProcData, "proc_layer");
    TTreeReaderValue<std::vector<int>> proc_time1(readerProcData, "proc_time1");
    TTreeReaderValue<std::vector<int>> proc_time2(readerProcData, "proc_time2");
    TTreeReaderValue<std::vector<bool>> has_valid_time_eta1(readerProcData, "has_valid_time_eta1");
    TTreeReaderValue<std::vector<bool>> has_valid_time_eta2(readerProcData, "has_valid_time_eta2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(readerTrackData, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(readerTrackData, "in_valid_track_eta2");

    long long event_count = 0;

    // Global counters
    long long total_hits_eta1 = 0;
    long long total_hits_eta2 = 0;
    long long total_hits_or = 0;
    long long total_hits_and = 0;

    // Eta/Layer-specific counters
    std::array<long long, 3> total_hits_per_layer_eta1 = {0, 0, 0};
    std::array<long long, 3> total_hits_per_layer_eta2 = {0, 0, 0};
    std::array<long long, 3> total_hits_per_layer_or = {0, 0, 0};
    std::array<long long, 3> total_hits_per_layer_and = {0, 0, 0};

    // Eta/Layer/Strip-specific counters
    long long hits_per_strip_eta1[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};
    long long hits_per_strip_eta2[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};
    long long hits_per_strip_or[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};
    long long hits_per_strip_and[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};
    while (readerProcData.Next() && readerTrackData.Next()) {
        event_count += 1;

        // Process Rates
        if (proc_time1.GetSetupStatus() == 0 && proc_time2.GetSetupStatus() == 0 &&
            has_valid_time_eta1.GetSetupStatus() == 0 && has_valid_time_eta2.GetSetupStatus() == 0)
        {
            const size_t n_hits = proc_time1->size();
            for (size_t i = 0; i < n_hits; ++i) {
                int layer = (*layers)[i];
                int strip = perFileHelpers::remapStrip((*strips)[i]);

                bool pass_track_cut = !in_valid_track_only || ((*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i]);
                if (!pass_track_cut) continue;

                if ((*has_valid_time_eta1)[i]) {
                    total_hits_eta1++;
                    total_hits_per_layer_eta1[layer]++;
                    hits_per_strip_eta1[layer][strip]++;
                }
                if ((*has_valid_time_eta2)[i]) {
                    total_hits_eta2++;
                    total_hits_per_layer_eta2[layer]++;
                    hits_per_strip_eta2[layer][strip]++;
                }
                if ((*has_valid_time_eta1)[i] || (*has_valid_time_eta2)[i]) {
                    total_hits_or++;
                    total_hits_per_layer_or[layer]++;
                    hits_per_strip_or[layer][strip]++;
                }
                if ((*has_valid_time_eta1)[i] && (*has_valid_time_eta2)[i]) {
                    total_hits_and++;
                    total_hits_per_layer_and[layer]++;
                    hits_per_strip_and[layer][strip]++;
                }
            }
        }
    }

    // Calculate final global rates
    auto calculateRate = [](long long hits, long long events, int n_strips, int n_layers)
        -> std::pair<double, double> {

        if (events <= 0 || n_strips <= 0 || n_layers <= 0) {
            return {0.0, 0.0};
        }

        const double total_exposure = events * TRIGGER_TIME_WINDOW * n_layers * n_strips * STRIP_WIDTH_CM * DETECTOR_LENGTH_CM;

        const double rate = static_cast<double>(hits) / total_exposure;
        const double rate_error = (hits > 0) ? (std::sqrt(static_cast<double>(hits)) / total_exposure) : 0.0;

        return {rate, rate_error};
    };

    // Calculate detector-wide rates
    std::tie(rate_results.rate_eta1, rate_results.rate_eta1_error) =
            calculateRate(total_hits_eta1, event_count, STRIPS_PER_LAYER, LAYER_COUNT);
    std::tie(rate_results.rate_eta2, rate_results.rate_eta2_error) =
            calculateRate(total_hits_eta2, event_count, STRIPS_PER_LAYER, LAYER_COUNT);
    std::tie(rate_results.rate_or, rate_results.rate_or_error) =
            calculateRate(total_hits_or, event_count, STRIPS_PER_LAYER, LAYER_COUNT);
    std::tie(rate_results.rate_and, rate_results.rate_and_error) =
            calculateRate(total_hits_and, event_count, STRIPS_PER_LAYER, LAYER_COUNT);

    // Calculate final Eta/Layer-specific rates
    for (int layer = 0; layer < LAYER_COUNT; ++layer) {

        std::tie(rate_results.rate_layers_eta1[layer], rate_results.rate_layers_eta1_error[layer]) =
            calculateRate(total_hits_per_layer_eta1[layer], event_count, STRIPS_PER_LAYER, 1);

        std::tie(rate_results.rate_layers_eta2[layer], rate_results.rate_layers_eta2_error[layer]) =
            calculateRate(total_hits_per_layer_eta2[layer], event_count, STRIPS_PER_LAYER, 1);

        std::tie(rate_results.rate_layers_or[layer], rate_results.rate_layers_or_error[layer]) =
            calculateRate(total_hits_per_layer_or[layer], event_count, STRIPS_PER_LAYER, 1);

        std::tie(rate_results.rate_layers_and[layer], rate_results.rate_layers_and_error[layer]) =
            calculateRate(total_hits_per_layer_and[layer], event_count, STRIPS_PER_LAYER, 1);

        // Calculate final Strip-specific rates
        for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {
            std::tie(rate_results.rate_strips_eta1[layer][strip],
                     rate_results.rate_strips_eta1_error[layer][strip]) =
                calculateRate(hits_per_strip_eta1[layer][strip], event_count, 1, 1);

            std::tie(rate_results.rate_strips_eta2[layer][strip],
                     rate_results.rate_strips_eta2_error[layer][strip]) =
                calculateRate(hits_per_strip_eta2[layer][strip], event_count, 1, 1);

            std::tie(rate_results.rate_strips_or[layer][strip],
                     rate_results.rate_strips_or_error[layer][strip]) =
                calculateRate(hits_per_strip_or[layer][strip], event_count, 1, 1);

            std::tie(rate_results.rate_strips_and[layer][strip],
                     rate_results.rate_strips_and_error[layer][strip]) =
                calculateRate(hits_per_strip_and[layer][strip], event_count, 1, 1);
        }
    }
}

void getDeadStrips(DeadStrips& dead_strips) {
    for (auto [layer, strip, side] : DEAD_STRIPS) {
        if (side == 1) {
            dead_strips.dead_strips_eta1.push_back({layer, strip});
        } else {
            dead_strips.dead_strips_eta2.push_back({layer, strip});
        }
    }
}

void getAverageToT(TFile* input_file, ToTResults& tot_results, bool in_beam_only,
    bool in_valid_track_only, DeadStrips& dead_strips, const RoI& region_of_interest)
{
    TTree* processed_tree = input_file->Get<TTree>("ProcessedData");
    TTree* track_tree = input_file->Get<TTree>("TrackReconstruction");

    if (!processed_tree || !track_tree) {
        std::cerr << "Error: Required trees not found in file " << input_file->GetName() << "\n";
        std::exit(EXIT_FAILURE);
    }

    TTreeReader reader_processed(processed_tree);
    TTreeReader reader_track(track_tree);

    TTreeReaderValue<std::vector<int>>  tot1(reader_processed, "proc_tot1");
    TTreeReaderValue<std::vector<int>>  tot2(reader_processed, "proc_tot2");
    TTreeReaderValue<std::vector<int>>  strips(reader_processed, "proc_strip");
    TTreeReaderValue<std::vector<int>>  layers(reader_processed, "proc_layer");
    TTreeReaderValue<std::vector<int>>  dts(reader_processed, "proc_dt_time1_time2");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(reader_track, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(reader_track, "in_valid_track_eta2");

    // Local accumulators
    Accumulator eta1[LAYER_COUNT][STRIPS_PER_LAYER];
    Accumulator eta2[LAYER_COUNT][STRIPS_PER_LAYER];

    while (reader_processed.Next() && reader_track.Next()) {
        if (tot1.GetSetupStatus() != 0 || in_valid_track_eta1.GetSetupStatus() != 0 || dts.GetSetupStatus() != 0) continue;

        const size_t n_hits = std::min({tot1->size(), tot2->size(), strips->size(), layers->size(),
            dts->size(), in_valid_track_eta1->size(), in_valid_track_eta2->size()});

        for (size_t i = 0; i < n_hits; ++i) {
            int layer = (*layers)[i];
            int strip = perFileHelpers::remapStrip((*strips)[i]);
            int dt = (*dts)[i];

            if (layer < 0 || layer >= LAYER_COUNT || strip < 0 || strip >= STRIPS_PER_LAYER) continue;

            // Track Cut
            bool pass_track_cut = !in_valid_track_only || ((*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i]);

            // Beam Cut
            bool pass_beam_cut = true;
            if (in_beam_only) {
                pass_beam_cut = std::find(region_of_interest.beam_region.begin(),
                                          region_of_interest.beam_region.end(),
                                          std::make_pair(strip, dt)) != region_of_interest.beam_region.end();
            }

            // Accumulate if both beam and track cuts pass
            if (pass_track_cut && pass_beam_cut) {
                if ((*tot1)[i] > 0) {
                    eta1[layer][strip].sum += (*tot1)[i];
                    eta1[layer][strip].hits++;
                }
                if ((*tot2)[i] > 0) {
                    eta2[layer][strip].sum += (*tot2)[i];
                    eta2[layer][strip].hits++;
                }
            }
        }
    }

    auto isDead = [](int l, int s, const std::vector<std::pair<int, int>>& dead_list) {
        return std::find(dead_list.begin(), dead_list.end(), std::make_pair(l, s)) != dead_list.end();
    };

    auto assignAverageAndError = [](const Accumulator& acc, bool is_dead, double& avg, ErrorRange& err) {
        if (is_dead || acc.hits == 0) {
            avg = std::numeric_limits<double>::quiet_NaN();
            err = ErrorRange(std::numeric_limits<double>::quiet_NaN());
        } else {
            avg = static_cast<double>(acc.sum) / acc.hits;
            err = ErrorRange(std::sqrt(avg) / acc.hits);
        }
    };

    for (int layer = 0; layer < LAYER_COUNT; ++layer) {
        // Accumulators for layer-wide totals
        Accumulator layer_acc_eta1 = {0, 0};
        Accumulator layer_acc_eta2 = {0, 0};

        for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {
            bool dead_eta1 = isDead(layer, strip, dead_strips.dead_strips_eta1);
            bool dead_eta2 = isDead(layer, strip, dead_strips.dead_strips_eta2);

            // Assign individual strip averages & errors
            assignAverageAndError(eta1[layer][strip],
                                  dead_eta1,
                                  tot_results.avg_tot_eta1[layer][strip],
                                  tot_results.avg_tot_eta1_error[layer][strip]);

            assignAverageAndError(eta2[layer][strip],
                                  dead_eta2,
                                  tot_results.avg_tot_eta2[layer][strip],
                                  tot_results.avg_tot_eta2_error[layer][strip]);

            // Aggregate active strip metrics for layer-wide calculation
            if (!dead_eta1) {
                layer_acc_eta1.sum += eta1[layer][strip].sum;
                layer_acc_eta1.hits += eta1[layer][strip].hits;
            }
            if (!dead_eta2) {
                layer_acc_eta2.sum += eta2[layer][strip].sum;
                layer_acc_eta2.hits += eta2[layer][strip].hits;
            }
        }

        // Assign layer-wide averages & errors
        assignAverageAndError(layer_acc_eta1,
                              false, // Layer itself isn't "dead", check hits instead via assignAverageAndError logic
                              tot_results.avg_tot_eta1_layer[layer],
                              tot_results.avg_tot_eta1_layer_error[layer]);

        assignAverageAndError(layer_acc_eta2,
                              false,
                              tot_results.avg_tot_eta2_layer[layer],
                              tot_results.avg_tot_eta2_layer_error[layer]);
    }
}

void getAverageMultiplicity(TFile* input_file, MultiplicityResults& mult_results,
    bool in_valid_track_only, DeadStrips& dead_strips)
{
    TTree* input_tree = input_file->Get<TTree>("InputData");
    TTree* processed_tree = input_file->Get<TTree>("ProcessedData");
    TTree* track_tree = input_file->Get<TTree>("TrackReconstruction");

    if (!input_tree || !processed_tree || !track_tree) {
        std::cerr << "Error: Required trees not found in file " << input_file->GetName() << "\n";
        std::exit(EXIT_FAILURE);
    }

    TTreeReader reader_input(input_tree);
    TTreeReader reader_processed(processed_tree);
    TTreeReader reader_track(track_tree);

    TTreeReaderValue<std::vector<int>>  raw_time1(reader_input, "hit_raw_time1");
    TTreeReaderValue<std::vector<int>>  raw_time2(reader_input, "hit_raw_time2");
    TTreeReaderValue<std::vector<int>>  rise(reader_input, "hit_rise");
    TTreeReaderValue<std::vector<int>>  strips(reader_processed, "proc_strip");
    TTreeReaderValue<std::vector<int>>  layers(reader_processed, "proc_layer");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta1(reader_track, "in_valid_track_eta1");
    TTreeReaderValue<std::vector<bool>> in_valid_track_eta2(reader_track, "in_valid_track_eta2");

    // Local counters
    int n_events = 0;
    int eta1_per_layer_strip[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};
    int eta2_per_layer_strip[LAYER_COUNT][STRIPS_PER_LAYER] = {{0}};

    // Track how many events actually had at least 1 hit on a specific strip
    int active_events_eta1[LAYER_COUNT][STRIPS_PER_LAYER] = {0};
    int active_events_eta2[LAYER_COUNT][STRIPS_PER_LAYER] = {0};

    while (reader_input.Next() && reader_processed.Next() && reader_track.Next()) {
        if (raw_time1.GetSetupStatus() < 0 || raw_time2.GetSetupStatus() < 0 || rise.GetSetupStatus() < 0 ||
            strips.GetSetupStatus() < 0 || layers.GetSetupStatus() < 0 ||
            in_valid_track_eta1.GetSetupStatus() < 0 || in_valid_track_eta2.GetSetupStatus() < 0)
        {
            continue;
        }

        const size_t n_hits = std::min({
            raw_time1->size(), raw_time2->size(), rise->size(),
            strips->size(), layers->size(),
            in_valid_track_eta1->size(), in_valid_track_eta2->size()
        });

        // Temporary flags to see if a strip fired in this event
        bool hit_in_event_eta1[LAYER_COUNT][STRIPS_PER_LAYER] = {false};
        bool hit_in_event_eta2[LAYER_COUNT][STRIPS_PER_LAYER] = {false};

        for (size_t i = 0; i < n_hits; ++i) {
            if ((*rise)[i] == 0) continue; // Skip falling edge hits

            int layer = (*layers)[i];
            int strip = perFileHelpers::remapStrip((*strips)[i]);

            if (layer < 0 || layer >= LAYER_COUNT || strip < 0 || strip >= STRIPS_PER_LAYER) continue;

            // Accumulate for all hits or only for hits that are part of valid tracks based on the flag
            if (!in_valid_track_only) {
                if ((*raw_time1)[i] != 0) {
                    eta1_per_layer_strip[layer][strip]++;
                    hit_in_event_eta1[layer][strip] = true;
                }
                if ((*raw_time2)[i] != 0) {
                    eta2_per_layer_strip[layer][strip]++;
                    hit_in_event_eta2[layer][strip] = true;
                }
            } else {
                if ((*in_valid_track_eta1)[i] || (*in_valid_track_eta2)[i]) {
                    if ((*raw_time1)[i] != 0) {
                        eta1_per_layer_strip[layer][strip]++;
                        hit_in_event_eta1[layer][strip] = true;
                    }
                    if ((*raw_time2)[i] != 0) {
                        eta2_per_layer_strip[layer][strip]++;
                        hit_in_event_eta2[layer][strip] = true;
                    }
                }
            }
        }
        for (int layer = 0; layer < LAYER_COUNT; ++layer) {
            for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {
                if (hit_in_event_eta1[layer][strip]) active_events_eta1[layer][strip]++;
                if (hit_in_event_eta2[layer][strip]) active_events_eta2[layer][strip]++;
            }
        }
        n_events++;
    }

    auto isDead = [](int l, int s, const std::vector<std::pair<int, int>>& dead_list) {
        return std::find(dead_list.begin(), dead_list.end(), std::make_pair(l, s)) != dead_list.end();
    };

    auto assignAverageAndError = [](int count, int active_events, bool is_dead, double& avg, ErrorRange& err) {
    if (is_dead || active_events == 0) {
        avg = std::numeric_limits<double>::quiet_NaN();
        err = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    } else {
        avg = static_cast<double>(count) / active_events;
        double raw_err = std::sqrt(count) / active_events;

        double err_low = (avg - raw_err < 1.0) ? (avg - 1.0) : raw_err;
        double err_high = raw_err;

        err = {err_low, err_high};
    }
};

    for (int layer = 0; layer < LAYER_COUNT; ++layer) {
        for (int strip = 0; strip < STRIPS_PER_LAYER; ++strip) {

            bool dead_eta1 = isDead(layer, strip, dead_strips.dead_strips_eta1);
            bool dead_eta2 = isDead(layer, strip, dead_strips.dead_strips_eta2);

            assignAverageAndError(eta1_per_layer_strip[layer][strip],
                                  active_events_eta1[layer][strip],
                                  dead_eta1,
                                  mult_results.avg_multiplicity_eta1[layer][strip],
                                  mult_results.avg_multiplicity_eta1_error[layer][strip]);

            assignAverageAndError(eta2_per_layer_strip[layer][strip],
                                  active_events_eta2[layer][strip],
                                  dead_eta2,
                                  mult_results.avg_multiplicity_eta2[layer][strip],
                                  mult_results.avg_multiplicity_eta2_error[layer][strip]);
        }
    }
}

void processToF(TFile* input_file, ToFResults& tof_results,
    TimeResolutionResults& time_resolution_results, DeadStrips& dead_strips)
{
    TTree* track_tree = input_file->Get<TTree>("TrackReconstruction");

    if (!track_tree) {
        std::cerr << "Error: Required trees not found in file " << input_file->GetName() << "\n";
        std::exit(EXIT_FAILURE);
    }

    TTreeReader reader_track(track_tree);

    // Global ToF Readers
    std::vector<std::unique_ptr<TTreeReaderValue<std::vector<int>>>> tof1_readers;
    std::vector<std::unique_ptr<TTreeReaderValue<std::vector<int>>>> tof2_readers;

    // Strip Readers
    std::vector<std::unique_ptr<TTreeReaderValue<std::vector<int>>>> stripF1_readers, stripS1_readers;
    std::vector<std::unique_ptr<TTreeReaderValue<std::vector<int>>>> stripF2_readers, stripS2_readers;

    // Data accumulation containers for percentile calculation
    std::array<std::vector<double>, LAYER_PAIR_COUNT> raw_global_eta1;
    std::array<std::vector<double>, LAYER_PAIR_COUNT> raw_global_eta2;

    std::array<std::array<std::vector<double>, STRIPS_PER_LAYER>, LAYER_PAIR_COUNT> raw_strip_eta1;
    std::array<std::array<std::vector<double>, STRIPS_PER_LAYER>, LAYER_PAIR_COUNT> raw_strip_eta2;

    for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {
        std::string base = "track_time_of_flight_layer_" + LAYER_PAIR_SUFFIXES[i];

        // Link Eta 1 branches
        tof1_readers.push_back(std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_track, (base + "_eta1").c_str()));
        stripF1_readers.push_back(std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_track, (base + "_stripFirst_eta1").c_str()));
        stripS1_readers.push_back(std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_track, (base + "_stripSecond_eta1").c_str()));

        // Link Eta 2 branches
        tof2_readers.push_back(std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_track, (base + "_eta2").c_str()));
        stripF2_readers.push_back(std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_track, (base + "_stripFirst_eta2").c_str()));
        stripS2_readers.push_back(std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_track, (base + "_stripSecond_eta2").c_str()));
    }

    // Pass 1: Accumulate raw values into memory (Decoupled branch verification)
    while (reader_track.Next()) {
        for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {

            // Process Eta 1 Global
            if (tof1_readers[i]->GetSetupStatus() == 0) {
                const auto& tofs = **tof1_readers[i];
                for (size_t k = 0; k < tofs.size(); ++k) {
                    double t_ns = TimeUtils::ticksToTime(tofs[k]);

                    raw_global_eta1[i].push_back(t_ns);
                    tof_results.time_of_flight_eta1[i].push_back(tofs[k]);

                    // Process Eta 1 Strips if available
                    if (stripF1_readers[i]->GetSetupStatus() == 0 &&
                        stripS1_readers[i]->GetSetupStatus() == 0)
                    {
                        const auto& sF = **stripF1_readers[i];
                        const auto& sS = **stripS1_readers[i];
                        if (k < sF.size() && k < sS.size()) {
                            int mapped_sF = perFileHelpers::remapStrip(sF[k]);
                            int mapped_sS = perFileHelpers::remapStrip(sS[k]);

                            if (mapped_sF == mapped_sS && mapped_sF >= 0 && mapped_sF < STRIPS_PER_LAYER) {
                                raw_strip_eta1[i][mapped_sF].push_back(t_ns);
                            }
                        }
                    }
                }
            }

            // Process Eta 2 Global
            if (tof2_readers[i]->GetSetupStatus() == 0) {
                const auto& tofs = **tof2_readers[i];
                for (size_t k = 0; k < tofs.size(); ++k) {
                    double t_ns = TimeUtils::ticksToTime(tofs[k]);

                    raw_global_eta2[i].push_back(t_ns);
                    tof_results.time_of_flight_eta2[i].push_back(tofs[k]);

                    // Process Eta 2 Strips if available
                    if (stripF2_readers[i]->GetSetupStatus() == 0 &&
                        stripS2_readers[i]->GetSetupStatus() == 0)
                    {
                        const auto& sF = **stripF2_readers[i];
                        const auto& sS = **stripS2_readers[i];
                        if (k < sF.size() && k < sS.size()) {
                            int mapped_sF = perFileHelpers::remapStrip(sF[k]);
                            int mapped_sS = perFileHelpers::remapStrip(sS[k]);

                            if (mapped_sF == mapped_sS && mapped_sF >= 0 && mapped_sF < STRIPS_PER_LAYER) {
                                raw_strip_eta2[i][mapped_sF].push_back(t_ns);
                            }
                        }
                    }
                }
            }
        }
    }

    // Reusable Lambda to perform the Two-Pass Gaussian Fit
    auto fitAndExtract = [](TH1D* hist, double& mean_out, ErrorRange& mean_err_out,
                            double& res_out, ErrorRange& res_err_out) {

        mean_out = std::numeric_limits<double>::quiet_NaN();
        mean_err_out = ErrorRange{std::numeric_limits<double>::quiet_NaN()};
        res_out = std::numeric_limits<double>::quiet_NaN();
        res_err_out = ErrorRange{std::numeric_limits<double>::quiet_NaN()};

        if (!hist || hist->Integral() < 20) return;

        // Pass 1: Global fit to find the general location of the peak
        TFitResultPtr r1 = hist->Fit("gaus", "Q0S");

        if (r1.Get() != nullptr && r1->IsValid() && static_cast<int>(r1) == 0) {
            double p_mean  = r1->Parameter(1);
            double p_sigma = r1->Parameter(2);
            double p_sigma_err = r1->Error(2);

            if (p_sigma <= 0.0) return;

            // Pass 2: Core fit restricted to ignore tails
            double fit_width = std::max(1.5 * p_sigma, 2.0 * TIME_TICK_NS);
            TFitResultPtr r2 = hist->Fit("gaus", "Q0S", "", p_mean - fit_width, p_mean + fit_width);

            if (r2.Get() != nullptr && r2->IsValid() && static_cast<int>(r2) == 0) {
                mean_out = r2->Parameter(1);
                mean_err_out = ErrorRange{r2->Error(1)};
                p_sigma = r2->Parameter(2);
                p_sigma_err = r2->Error(2);
            } else {
                mean_out = p_mean;
                mean_err_out = ErrorRange{r1->Error(1)};
            }

            const double TDC_VARIANCE = (TIME_TICK_NS * TIME_TICK_NS) / 12.0;

            double raw_res_sq = (p_sigma * p_sigma) / 2.0;
            double intrinsic_res_sq = raw_res_sq - TDC_VARIANCE;

            if (intrinsic_res_sq > 0.0) {
                res_out = std::sqrt(intrinsic_res_sq);
                double propagated_err = (p_sigma / (2.0 * res_out)) * p_sigma_err;

                double err_low = (res_out - propagated_err < 0.0) ? res_out : propagated_err;
                double err_high = propagated_err;

                res_err_out = ErrorRange{err_low, err_high};
            }
        }
    };

    // Pass 2: Dynamically construct histograms using percentile ranges & run fits
    double bin_width = TIME_TICK_NS;
    for (int i = 0; i < LAYER_PAIR_COUNT; ++i) {

        // Global Eta1/Eta2 Percentile Dynamic Boundaries (95% central interval)
        auto [g_xmin1, g_xmax1] = perFileHelpers::getPercentileRange(raw_global_eta1[i], 0.025, 0.975);
        auto [g_xmin2, g_xmax2] = perFileHelpers::getPercentileRange(raw_global_eta2[i], 0.025, 0.975);

        int n_bins1 = std::max(1, static_cast<int>(std::round((g_xmax1 - g_xmin1) / bin_width)));
        int n_bins2 = std::max(1, static_cast<int>(std::round((g_xmax2 - g_xmin2) / bin_width)));

        g_xmax1 = g_xmin1 + n_bins1 * bin_width;
        g_xmax2 = g_xmin2 + n_bins2 * bin_width;

        TH1D* h_tof1 = new TH1D(Form("h_tof1_%s", LAYER_PAIR_SUFFIXES[i].c_str()), Form("ToF Eta1 Layer %s", LAYER_PAIR_SUFFIXES[i].c_str()), n_bins1, g_xmin1, g_xmax1);
        TH1D* h_tof2 = new TH1D(Form("h_tof2_%s", LAYER_PAIR_SUFFIXES[i].c_str()), Form("ToF Eta2 Layer %s", LAYER_PAIR_SUFFIXES[i].c_str()), n_bins2, g_xmin2, g_xmax2);

        for (double val : raw_global_eta1[i]) h_tof1->Fill(val);
        for (double val : raw_global_eta2[i]) h_tof2->Fill(val);

        fitAndExtract(h_tof1,
                      tof_results.avg_time_of_flight_eta1[i], tof_results.avg_time_of_flight_eta1_error[i],
                      time_resolution_results.time_resolution_global_eta1[i], time_resolution_results.time_resolution_global_eta1_error[i]);

        fitAndExtract(h_tof2,
                      tof_results.avg_time_of_flight_eta2[i], tof_results.avg_time_of_flight_eta2_error[i],
                      time_resolution_results.time_resolution_global_eta2[i], time_resolution_results.time_resolution_global_eta2_error[i]);

        delete h_tof1;
        delete h_tof2;

        // Per-Strip Percentile Dynamic Boundaries & Fits
        for (int s = 0; s < STRIPS_PER_LAYER; ++s) {
            auto [s_xmin1, s_xmax1] = perFileHelpers::getPercentileRange(raw_strip_eta1[i][s], 0.025, 0.975);
            auto [s_xmin2, s_xmax2] = perFileHelpers::getPercentileRange(raw_strip_eta2[i][s], 0.025, 0.975);

            int n_bins1 = std::max(1, static_cast<int>(std::round((s_xmax1 - s_xmin1) / bin_width)));
            int n_bins2 = std::max(1, static_cast<int>(std::round((s_xmax2 - s_xmin2) / bin_width)));

            s_xmax1 = s_xmin1 + n_bins1 * bin_width;
            s_xmax2 = s_xmin2 + n_bins2 * bin_width;

            TH1D* h_strip1 = new TH1D(Form("h_tof1_%s_s%d", LAYER_PAIR_SUFFIXES[i].c_str(), s), "ToF", n_bins1, s_xmin1, s_xmax1);
            TH1D* h_strip2 = new TH1D(Form("h_tof2_%s_s%d", LAYER_PAIR_SUFFIXES[i].c_str(), s), "ToF", n_bins2, s_xmin2, s_xmax2);

            for (double val : raw_strip_eta1[i][s]) h_strip1->Fill(val);
            for (double val : raw_strip_eta2[i][s]) h_strip2->Fill(val);

            double dummy_mean_eta1 = 0.0, dummy_mean_eta2 = 0.0;
            ErrorRange dummy_mean_err_eta1 = {0.0, 0.0}, dummy_mean_err_eta2 = {0.0, 0.0};

            fitAndExtract(h_strip1,
                          dummy_mean_eta1, dummy_mean_err_eta1,
                          time_resolution_results.time_resolution_strip_eta1[i][s],
                          time_resolution_results.time_resolution_strip_eta1_error[i][s]);

            fitAndExtract(h_strip2,
                          dummy_mean_eta2, dummy_mean_err_eta2,
                          time_resolution_results.time_resolution_strip_eta2[i][s],
                          time_resolution_results.time_resolution_strip_eta2_error[i][s]);

            delete h_strip1;
            delete h_strip2;
        }
    }
}

}