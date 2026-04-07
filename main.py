# =============================================================================
# Interactive Bokeh app for visualizing the response of an SDOF system
# under Gaussian random excitation with different covariance structures.
#
# Run locally with:
#     bokeh serve --show DIRECTORY OF THIS FILE
# =============================================================================


# =============================================================================
# Imports
# =============================================================================

# math and scientific computing
import inspect
import statistics
import sys
from os.path import abspath, dirname, join

import numpy as np
from scipy.fft import rfftfreq
from scipy.stats import norm

# bokeh
from bokeh.io import curdoc
from bokeh.layouts import row

# latex support
currentdir = dirname(abspath(inspect.getfile(inspect.currentframe())))
parentdir = join(dirname(currentdir), "shared/")
sys.path.insert(0, parentdir)

# project modules
from calculations import (
    Properties,
    harmonic_transfer_function,
    prob_of_failure_time,
    prob_of_surv,
    rate_of_upcrossing,
    unit_impulse_function,
)
from gaussian_process import (
    Gaussian_distribution_cdf,
    Gaussian_distribution_pdf,
    cos_cov,
    dirac_delta_cov,
    gaussian_function,
    matérn_class,
    squared_exponential_cov,
)
from rfourier_utils import to_frequency_domain, to_time_domain
from ui import create_widgets, create_plots, layout


# =============================================================================
# Global constants and initial state
# =============================================================================

sample_rate = 128
observation_duration = 8
time_step = 1 / sample_rate
number_of_samples = int(observation_duration * sample_rate)

frequency_bins = rfftfreq(number_of_samples, time_step)
time_values = np.linspace(
    -observation_duration / 2,
    observation_duration / 2,
    number_of_samples,
    endpoint=False,
)

system_params = Properties()
system_params.stiffness = 10 * pow(20 * np.pi, 2)
system_params.mass = 10
system_params.zeta = 0.05


class Status:  # save the current values of excitation and response
    def __init__(self, n):
        self.excitation = None
        self.additional_realizations = []
        self.standard_deviation = 0
        self.cov = np.zeros(len(n))
        self.excitation_spectral_density = []
        self.response_time_series = None


status_values = Status(time_values)
# standard_deviation = 0
# cov = np.zeros(len(time_values))
# additional_realizations = []
# excitation_spectral_density = []
active_view = [1]
active_covariance_category = []


# =============================================================================
# Set up UI
# =============================================================================

widgets = create_widgets(system_params)
plot_components = create_plots(status_values.standard_deviation)

figures = plot_components["figures"]
sources = plot_components["sources"]

explain, visualization, final = layout(widgets, figures)

# =============================================================================
# Helper functions
# =============================================================================

def get_current_covariance_matrix():
    correlation = widgets["correlation_slider"].value
    sigma = widgets["sigma_slider"].value
    bandwidth = widgets["bandwidth_slider"].value
    freq = widgets["frequency_slider"].value

    if widgets["covariance_checkbox"].active == [0]:
        return cos_cov(time_values, freq, bandwidth, sigma)
    if widgets["covariance_checkbox"].active == [1]:
        return squared_exponential_cov(time_values, sigma, correlation)
    if widgets["covariance_checkbox"].active == [2]:
        return matérn_class(time_values, correlation, sigma)
    if widgets["covariance_checkbox"].active == [3]:
        return dirac_delta_cov(sigma, time_values)

    return np.zeros((len(time_values), len(time_values)))

def clear_source_xy(source):
    source.data = {"x": time_values, "y": []}

def excitation_real(time_values):
    """Generate one realization of the currently selected Gaussian excitation process."""
    mean = widgets["mean_slider"].value
    covariance_matrix = get_current_covariance_matrix()
    return gaussian_function(mean * np.ones(len(time_values)), covariance_matrix)


# =============================================================================
# Callback functions
# =============================================================================

def update_statistic_checkbox(attrname, old, new):
    """Toggle between visualization mode and theory mode."""
    global active_view

    if widgets["visu_or_theory_checkbox"].active == [0, 1]:
        if active_view == [0]:
            widgets["visu_or_theory_checkbox"].active = [1]
        else:
            widgets["visu_or_theory_checkbox"].active = [0]

    if widgets["visu_or_theory_checkbox"].active == [1]:
        widgets["covariance_checkbox"].visible = False
        visualization.visible = False
        explain.visible = True
        widgets["covariance_checkbox"].active = []

    if widgets["visu_or_theory_checkbox"].active == [0]:
        widgets["covariance_checkbox"].visible = True
        visualization.visible = True
        explain.visible = False

    active_view = widgets["visu_or_theory_checkbox"].active


def update_covariance_checkbox(attrname, old, new):
    """Ensure that only one covariance category is active at a time."""
    global active_covariance_category

    if widgets["covariance_checkbox"].active == [0, 1]:
        if active_covariance_category == [0]:
            widgets["covariance_checkbox"].active = [1]
        else:
            widgets["covariance_checkbox"].active = [0]

    if widgets["covariance_checkbox"].active == [1, 2]:
        if active_covariance_category == [1]:
            widgets["covariance_checkbox"].active = [2]
        else:
            widgets["covariance_checkbox"].active = [1]

    if widgets["covariance_checkbox"].active == [2, 3]:
        if active_covariance_category == [2]:
            widgets["covariance_checkbox"].active = [3]
        else:
            widgets["covariance_checkbox"].active = [2]

    if widgets["covariance_checkbox"].active == [0, 3]:
        if active_covariance_category == [0]:
            widgets["covariance_checkbox"].active = [3]
        else:
            widgets["covariance_checkbox"].active = [0]

    if widgets["covariance_checkbox"].active == [1, 3]:
        if active_covariance_category == [1]:
            widgets["covariance_checkbox"].active = [3]
        else:
            widgets["covariance_checkbox"].active = [1]

    if widgets["covariance_checkbox"].active == [0, 2]:
        if active_covariance_category == [0]:
            widgets["covariance_checkbox"].active = [2]
        else:
            widgets["covariance_checkbox"].active = [0]

    update_cov_slider()
    active_covariance_category = widgets["covariance_checkbox"].active


def update_cov_slider():
    """Update visible sliders and reset response plots when covariance type changes."""
    for source in [
        sources["ds_response_time"],
        sources["ds_responsef_freq"],
        sources["ds_response_t_mean"],
        sources["ds_spectral_density_res"],
        sources["ds_covariance_res"],
        sources["ds_ex_new1"],
        sources["ds_ex_new2"],
        sources["ds_ex_new3"],
        sources["res_new_1"],
        sources["res_new_2"],
        sources["res_new_3"],
    ]:
        clear_source_xy(source)

    sources["ds_response_varea1"].data = {"x": time_values, "y1": [], "y2": []}
    sources["ds_response_varea2"].data = {"x": time_values, "y1": [], "y2": []}

    status_values.additional_realizations = []

    if widgets["covariance_checkbox"].active == [0]:
        widgets["sigma_slider"].visible = True
        widgets["mean_slider"].visible = True
        widgets["correlation_slider"].visible = False
        widgets["frequency_slider"].visible = True
        widgets["bandwidth_slider"].visible = True

    elif widgets["covariance_checkbox"].active == [1] or widgets["covariance_checkbox"].active == [2]:
        widgets["sigma_slider"].visible = True
        widgets["mean_slider"].visible = True
        widgets["correlation_slider"].visible = True
        widgets["frequency_slider"].visible = False
        widgets["bandwidth_slider"].visible = False

    elif widgets["covariance_checkbox"].active == [3]:
        widgets["sigma_slider"].visible = True
        widgets["mean_slider"].visible = True
        widgets["correlation_slider"].visible = False
        widgets["frequency_slider"].visible = False
        widgets["bandwidth_slider"].visible = False

    update_excitation_cov(None, None, None)


def update_impulse(attrname, old, new):
    """Update the SDOF impulse response and transfer function plots."""
    system_params.stiffness = widgets["stiffness_slider"].value
    system_params.mass = widgets["mass_slider"].value
    system_params.zeta = widgets["damping_slider"].value

    widgets["frequency_title"].text = (
        f"Natural frequency: {round(system_params.w, 2)}rad/s "
        f"or {round(system_params.w / (2 * np.pi), 2)}Hz"
    )

    h_x = unit_impulse_function(time_values, system_params)
    H_x = to_frequency_domain(time_values, h_x)

    sources["ds_transfer_t"].data = {"x": time_values, "y": h_x}
    sources["ds_transfer_f_mag"].data = {"x": frequency_bins, "y": np.abs(H_x)}


def update_excitation_cov(attrname, old, new):
    """Refresh excitation time-domain and frequency-domain plots after parameter changes."""
    status_values.excitation = excitation_real(time_values)
    excitation_freq_spectrum = to_frequency_domain(time_values, status_values.excitation)

    mean = widgets["mean_slider"].value
    sigma = widgets["sigma_slider"].value

    sources["ds_excitation_t"].data = {"x": time_values, "y": status_values.excitation}
    sources["ds_excitation_f_real"].data = {
        "x": frequency_bins,
        "y": np.abs(excitation_freq_spectrum),
    }

    sources["ds_excitation_time_varea"].data = {
        "x": time_values,
        "y1": (mean - sigma) * np.ones(len(time_values)),
        "y2": (mean + sigma) * np.ones(len(time_values)),
    }

    sources["ds_excitation_time_3_varea"].data = {
        "x": time_values,
        "y1": (mean - 3 * sigma) * np.ones(len(time_values)),
        "y2": (mean + 3 * sigma) * np.ones(len(time_values)),
    }

    sources["ds_excitation_t_mean"].data = {
        "x": time_values,
        "y": mean * np.ones(len(time_values)),
    }

    covariance_distribution()


def covariance_distribution():
    """Compute and plot the excitation covariance function for the selected process."""
    cov_matrix = get_current_covariance_matrix()
    center_index = len(time_values) // 2
    status_values.cov = cov_matrix[center_index]

    sources["ds_covariance_exc"].data = {"x": time_values, "y": status_values.cov}
    spec_dens(status_values.cov)


def spec_dens(cov):
    """Compute and plot the spectral density from the current covariance function."""
    status_values.excitation_spectral_density = to_frequency_domain(time_values, cov)
    sources["ds_spectral_density_exc"].data = {
        "x": frequency_bins,
        "y": np.abs(status_values.excitation_spectral_density),
    }


def add_realization():
    """Add one additional excitation realization to the plot."""
    status_values.additional_realizations.append(excitation_real(time_values))

    if len(status_values.additional_realizations) == 1:
        sources["ds_ex_new1"].data = {"x": time_values, "y": status_values.additional_realizations[0]}

    if len(status_values.additional_realizations) == 2:
        sources["ds_ex_new2"].data = {"x": time_values, "y": status_values.additional_realizations[1]}

    if len(status_values.additional_realizations) == 3:
        sources["ds_ex_new3"].data = {"x": time_values, "y": status_values.additional_realizations[2]}

    if len(status_values.additional_realizations) > 3:
        status_values.additional_realizations = status_values.additional_realizations[:3]


def remove_realization():
    """Remove the most recently added excitation realization."""
    if len(status_values.additional_realizations) == 1:
        clear_source_xy(sources["ds_ex_new1"])
        status_values.additional_realizations.pop()

    if len(status_values.additional_realizations) == 2:
        clear_source_xy(sources["ds_ex_new2"])
        status_values.additional_realizations.pop()

    if len(status_values.additional_realizations) == 3:
        clear_source_xy(sources["ds_ex_new3"])
        status_values.additional_realizations.pop()


def update_response():
    """Update the response realizations and then compute response statistics."""
    clear_source_xy(sources["res_new_1"])
    clear_source_xy(sources["res_new_2"])
    clear_source_xy(sources["res_new_3"])
    update_response_covariance()


def update_response_covariance():
    """Compute the SDOF response from excitation and update response plots."""
    system_params.stiffness = widgets["stiffness_slider"].value
    system_params.mass = widgets["mass_slider"].value
    system_params.zeta = widgets["damping_slider"].value

    H_x = harmonic_transfer_function(frequency_bins, system_params)
    excitation_freq_spectrum = to_frequency_domain(time_values, status_values.excitation)
    response_freq_spectrum = H_x * excitation_freq_spectrum
    status_values.response_time_series = to_time_domain(frequency_bins, response_freq_spectrum)

    sources["ds_response_time"].data = {"x": time_values, "y": status_values.response_time_series.real}
    sources["ds_responsef_freq"].data = {
        "x": frequency_bins,
        "y": np.abs(response_freq_spectrum),
    }

    if len(status_values.additional_realizations) == 1:
        excitation_freq_spectrum = to_frequency_domain(time_values, status_values.additional_realizations[0])
        response_freq_spectrum = H_x * excitation_freq_spectrum
        response_time_series = to_time_domain(frequency_bins, response_freq_spectrum)
        sources["res_new_1"].data = {"x": time_values, "y": response_time_series.real}

    if len(status_values.additional_realizations) == 2:
        excitation_freq_spectrum_1 = to_frequency_domain(time_values, status_values.additional_realizations[0])
        response_freq_spectrum_1 = H_x * excitation_freq_spectrum_1
        response_time_series_1 = to_time_domain(frequency_bins, response_freq_spectrum_1)

        excitation_freq_spectrum_2 = to_frequency_domain(time_values, status_values.additional_realizations[1])
        response_freq_spectrum_2 = H_x * excitation_freq_spectrum_2
        response_time_series_2 = to_time_domain(frequency_bins, response_freq_spectrum_2)

        sources["res_new_1"].data = {"x": time_values, "y": response_time_series_1.real}
        sources["res_new_2"].data = {"x": time_values, "y": response_time_series_2.real}

    if len(status_values.additional_realizations) == 3:
        excitation_freq_spectrum_1 = to_frequency_domain(time_values, status_values.additional_realizations[0])
        response_freq_spectrum_1 = H_x * excitation_freq_spectrum_1
        response_time_series_1 = to_time_domain(frequency_bins, response_freq_spectrum_1)

        excitation_freq_spectrum_2 = to_frequency_domain(time_values, status_values.additional_realizations[1])
        response_freq_spectrum_2 = H_x * excitation_freq_spectrum_2
        response_time_series_2 = to_time_domain(frequency_bins, response_freq_spectrum_2)

        excitation_freq_spectrum_3 = to_frequency_domain(time_values, status_values.additional_realizations[2])
        response_freq_spectrum_3 = H_x * excitation_freq_spectrum_3
        response_time_series_3 = to_time_domain(frequency_bins, response_freq_spectrum_3)

        sources["res_new_1"].data = {"x": time_values, "y": response_time_series_1.real}
        sources["res_new_2"].data = {"x": time_values, "y": response_time_series_2.real}
        sources["res_new_3"].data = {"x": time_values, "y": response_time_series_3.real}

    calculate_stochastic(status_values.cov, status_values.excitation_spectral_density)


def calculate_stochastic(excitation_covariance, excitation_spectral_density):
    """Calculate response covariance, spectral density, and statistical envelopes."""
    mean = statistics.mean(status_values.response_time_series)

    sources["ds_response_t_mean"].data = {
        "x": time_values,
        "y": np.ones(len(time_values)) * mean,
    }

    H_x = harmonic_transfer_function(frequency_bins, system_params)

    if widgets["covariance_checkbox"].active == [3]:
        excitation_covariance = excitation_covariance / time_step
        excitation_spectral_density = to_frequency_domain(time_values, excitation_covariance)

    response_spectral_density = np.multiply(np.abs(H_x) ** 2, excitation_spectral_density)

    sources["ds_spectral_density_res"].data = {
        "x": frequency_bins,
        "y": np.abs(response_spectral_density),
    }

    response_covariance = to_time_domain(frequency_bins, response_spectral_density)
    sources["ds_covariance_res"].data = {"x": time_values, "y": response_covariance.real}

    center_index = len(time_values) // 2
    variance = response_covariance.real[center_index]
    status_values.standard_deviation = np.sqrt(variance)

    three_sigma = 3 * status_values.standard_deviation

    sources["ds_response_varea1"].data = {
        "x": time_values,
        "y1": (mean - status_values.standard_deviation) * np.ones(len(time_values)),
        "y2": (mean + status_values.standard_deviation) * np.ones(len(time_values)),
    }

    sources["ds_response_varea2"].data = {
        "x": time_values,
        "y1": (mean - three_sigma) * np.ones(len(time_values)),
        "y2": (mean + three_sigma) * np.ones(len(time_values)),
    }

    update_gaussian_distribution()
    update_failure_analysis(None, None, None)


def update_gaussian_distribution():
    """Plot the marginal Gaussian PDF and CDF of the current response."""
    mean = statistics.mean(status_values.response_time_series)

    u = np.linspace(
        -5 * status_values.standard_deviation + mean,
        5 * status_values.standard_deviation + mean,
        100,
    )
    gauss_pdf = Gaussian_distribution_pdf(u, mean, status_values.standard_deviation)
    gauss_cdf = Gaussian_distribution_cdf(u, mean, status_values.standard_deviation)

    x_sigma = np.linspace(
        -status_values.standard_deviation + mean,
        status_values.standard_deviation + mean,
        20,
    )
    y_sigma = gauss_pdf[40:60]
    y_sigma_bottom = np.zeros(len(y_sigma))

    x_three_sigma = np.linspace(
        -3 * status_values.standard_deviation + mean,
        3 * status_values.standard_deviation + mean,
        60,
    )
    y_three_sigma = gauss_pdf[20:80]
    y_three_sigma_bottom = np.zeros(len(y_three_sigma))

    sources["ds_pdf_gaussian_vbar1"].data = {
        "x": x_sigma,
        "top": y_sigma,
        "bottom": y_sigma_bottom,
    }

    sources["ds_pdf_gaussian_vbar2"].data = {
        "x": x_three_sigma,
        "top": y_three_sigma,
        "bottom": y_three_sigma_bottom,
    }

    sources["ds_pdf_gaussian"].data = {"x": u, "y": gauss_pdf}
    sources["ds_cdf_gaussian"].data = {"x": u, "y": gauss_cdf}


def update_failure_analysis(attr, old, new):
    """Compute failure and survival probabilities for the chosen threshold."""
    mean = statistics.mean(status_values.response_time_series)
    velocity = np.gradient(status_values.response_time_series, time_step)
    sigma_velocity = np.std(velocity)

    u = widgets["failure_slider"].value * status_values.standard_deviation
    v = 2 * rate_of_upcrossing(u, status_values.standard_deviation, sigma_velocity, 0)

    d = norm(loc=mean, scale=status_values.standard_deviation)
    L = d.cdf(np.abs(u)) - d.cdf(np.abs(-u))

    prob_surv = prob_of_surv(v, L, time_values)
    prob_time = prob_of_failure_time(v, time_values)

    center_index = len(time_values) // 2

    sources["ds_prob_of_survival"].data = {
        "x": time_values[center_index:],
        "y": prob_surv[center_index:],
    }

    sources["ds_prob_of_time"].data = {
        "x": time_values[center_index:],
        "y": prob_time[center_index:],
    }


# =============================================================================
# Widget callback bindings
# =============================================================================

widgets["visu_or_theory_checkbox"].on_change("active", update_statistic_checkbox)
widgets["covariance_checkbox"].on_change("active", update_covariance_checkbox)

widgets["damping_slider"].on_change("value_throttled", update_impulse)
widgets["stiffness_slider"].on_change("value_throttled", update_impulse)
widgets["mass_slider"].on_change("value_throttled", update_impulse)

widgets["sigma_slider"].on_change("value_throttled", update_excitation_cov)
widgets["mean_slider"].on_change("value_throttled", update_excitation_cov)
widgets["correlation_slider"].on_change("value_throttled", update_excitation_cov)
widgets["frequency_slider"].on_change("value_throttled", update_excitation_cov)
widgets["bandwidth_slider"].on_change("value_throttled", update_excitation_cov)

widgets["add_gaussian_button"].on_click(add_realization)
widgets["remove_gaussian_button"].on_click(remove_realization)

widgets["calculate_button_gauss"].on_click(update_response)

widgets["failure_slider"].on_change("value_throttled", update_failure_analysis)


# =============================================================================
# Initial app setup
# =============================================================================

update_impulse(None, None, None)
update_statistic_checkbox(None, None, None)
covariance_distribution()
spec_dens(status_values.cov)

# =============================================================================
# Bokeh document setup
# =============================================================================

curdoc().add_root(row(final, width=3000))
curdoc().title = "SDOF under Random Vibration"