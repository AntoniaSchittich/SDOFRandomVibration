# math and scientific computing
import numpy as np

# bokeh
from bokeh.models import (
    Button,
    CheckboxButtonGroup,
    Div,
    Slider,
    PrintfTickFormatter,
)

from bokeh.plotting import figure
from bokeh.layouts import column, row


from latex_support import LatexDiv, LatexLabel

plot_width = 300
plot_height = plot_width * 2 // 3
line_width = 1
text_font_size = "10pt"

def style_figure(fig, scientific_y=False):
    fig.title.align = "center"
    fig.toolbar.logo = None
    fig.outline_line_width = line_width
    fig.outline_line_color = "Black"
    fig.title.text_font_size = text_font_size
    if scientific_y:
        fig.yaxis[0].formatter = PrintfTickFormatter(format="%4.0e")



def create_widgets(system_params):
    
    # checkbox groups
    visu_or_theory_checkbox = CheckboxButtonGroup(
        labels=["Process Category", "Theory"],
        active=[1],
        width=1250,
    )

    covariance_checkbox = CheckboxButtonGroup(
        labels=[
            "Narrowband Process",
            " Squared Exponential Cov",
            "Martérn Class",
            "White Noise",
        ],
        active=[],
        width=1250,
    )

    # buttons
    calculate_button_Gauss = Button(label="Calculate", button_type="primary", height=40)
    Add_Gaussian_Button = Button(label="Add Realization", button_type="primary", height=30)
    Remove_Gaussian_Button = Button(label="Remove Realization", button_type="primary", height=30)

    # sliders
    damping_slider = Slider(
        title="damping ratio \u03B6 [-]",
        value=system_params.zeta,
        start=0.01,
        end=1,
        step=0.01,
    )

    stiffness_slider = Slider(
        title="stiffness k [N/m]",
        value=system_params.stiffness,
        start=1000,
        end=100000,
        step=1000,
    )

    mass_slider = Slider(
        title="mass m [kg]",
        value=system_params.mass,
        start=1,
        end=25,
        step=1,
    )

    sigma_slider = Slider(
        title="standard deviation \u03C3 [-]",
        value=2,
        start=0.0,
        end=5.0,
        step=0.5,
        max_width=300,
    )

    mean_slider = Slider(
        title="mean value \u03BC [-]",
        value=0,
        start=-5,
        end=5.0,
        step=1,
    )

    correlation_slider = Slider(
        title="correlation lenght/  length scale [-]",
        value=0.1,
        start=0.001,
        end=1,
        step=0.001,
        max_width=300,
    )

    bandwidth_slider = Slider(
        title="bandwidth [-]",
        value=0,
        start=0,
        end=np.pi * 2,
        step=np.pi / 8,
        max_width=300,
    )

    frequency_slider = Slider(
        title="frequency f [Hz]",
        value=5,
        start=1,
        end=10,
        step=1,
        max_width=300,
    )

    failure_slider = Slider(
        title="critical level |u|[-]",
        value=3,
        start=1,
        end=6,
        step=0.5,
        min_width=600,
        max_width=600,
    )

    div_width = 300

    title_div = Div(
        text="<b>SDOF under Random Vibration</b>",
        style={"font-size": "30px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    excitation_div = Div(
        text="<b>Excitation</b>",
        style={"font-size": "20px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    para_div = Div(
        text="<b>Parameter Control</b>",
        style={"font-size": "20px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
    )

    Response_div = Div(
        text="<b>Response</b>",
        style={"font-size": "20px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    Gaussian_div = Div(
        text="<b>Marginal Gaussian Distribution</b>",
        style={"font-size": "20px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    Control_para_div = Div(
        text="<b>Control Parameters</b>",
        style={"font-size": "20px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    Failure_div = Div(
        text="<b>Failure Analysis</b>",
        style={"font-size": "20px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    system_para_div = Div(
        text="<b>SDOF System</b>",
        style={"font-size": "15px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    SDOF_system_div = Div(
        text="<b>SDOF System</b>",
        style={"font-size": "20px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    excitation_para = Div(
        text="<b>Excitation</b>",
        style={"font-size": "15px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    realization_div = Div(
        text="<b>Realization</b>",
        style={"font-size": "15px"},
        width=div_width,
        height=15,
        sizing_mode="stretch_width",
        align="center",
        height_policy="fit",
    )

    frequency_title = Div(
        text=f"Natural frequency: {round(system_params.w, 2)}rad/s or {round(system_params.w / (2 * np.pi), 2)}Hz",
        style={"font-size": "12px"},
    )

    theory_div = LatexDiv(
        text=
        " This application is supposed to visualize the response of a single degree of freedom (SDOF) system subjected to random vibration. "
        "The excitation thereby is a Gaussian random process $\\{X(t)\\}$ and thus solely dependent on the mean vector "
        "$\\bm{\\mu}_{X}$ and covariance matrix $\\mathbf{K}_{X X}$:"
        "$$p_{\\mathbf{X}}(\\mathbf{u})=\\frac{1}{(2\\pi)^{n / 2}\\left|\\mathbf{K}_{X X}\\right|^{1 / 2}} "
        "\\exp \\left(-\\frac{1}{2}\\left(\\mathbf{u}-\\bm{\\mu}_{X}\\right)^{T}\\mathbf{K}_{X X}^{-1}"
        "\\left(\\mathbf{u}-\\bm{\\mu}_{X}\\right)\\right)\\:.$$"
        "The response of the SDOF system is calculated with the harmonic transfer function $H_x(f)$ in the frequency domain "
        "with the following formular:"
        "$$X(f) = H_x(f)F(f)\\:,$$"
        "where $F(f)$ and $X(f)$ denote the excitation and response spectrum, respectively."
        "</p>"
        "Furthermore, all processes are second-order stationary. This indicates that the mean and the covarinace are only dependent "
        "on the time difference between the observations along the abscissa. Thus, the mean value is constant for all values of $t$:"
        "$$\\mu_{X}(t)=\\mu_{X}\\:.$$"
        "The reason for considering stationary vibrations is that we suppose the system has oscillated since $t = -\\infty$. "
        "Thus, we can assume a steady-state response. Next to a good approximation, the choice of a Gaussian random process enables "
        "a feasible response calculation. As it is known that the response of a system subjected to a Gaussian vibration will also "
        "be Gaussian distributed."
        "</p>"
        "On top of the application, the user can choose between different categories of random processes. The categories are defined "
        "by their covaraince and thus by the degree of correlations between the different random variables that define the Gaussian process. "
        "Possible choices and their respective stationary covariance functions are:"
        "</p>"
        "Narrowband Processes: $$G_{\\mathrm{XX,NB}}(\\tau) = \\sigma^2 \\cos{(f_c\\tau)} "
        "\\frac{\\sin{(b\\tau)}}{b\\tau}\\:,$$"
        "Processes with a Squared Exponential Covariance: $$G_{\\mathrm{XX,SE}}(\\tau)= "
        "\\sigma^2 \\exp \\left(-\\frac{\\tau^{2}}{2 \\ell^{2}}\\right)\\:,$$"
        "Processes with a Matérn Class Covariance :$$G_{\\mathrm{XX,MC}}(\\tau)= "
        "\\sigma^{2}\\frac{2^{1-\\nu}}{\\Gamma(\\nu)}\\left(\\frac{\\sqrt{2 \\nu} \\tau}{\\ell}\\right)^{\\nu} "
        "K_{\\nu}\\left(\\frac{\\sqrt{2 \\nu} \\tau}{\\ell}\\right) \\:,$$ and "
        "White Noise Processes: $$  G_{\\mathrm{XX,WN}}(\\tau) = 2\\pi S_0\\delta(\\tau)\\:,$$"
        "where \\tau describes the time difference between the observations, \\sigma the standard deviation, "
        "f_c the characteristic frequency, b the bandwidth of the spectrum, \\ell the correlation length, "
        "$\\Gamma(\\cdot)$ and $K_{\\nu}(\\cdot)$ represent the gamma function and the Bessel function, "
        "and S_0 the autospectral density level."
        "</p>"
        "After choosing a covariance the system and the excitation can be altered by the control parameters on the top of the page. "
        "To demonstrate that the shown vibration is only one possible realization, the user can add and delete other additional_realizations."
        "</p>"
        "At the bottom of the application, the marginal cumulative distribution function (CDF) and marginal probability density function (PDF) "
        "of each random variable of the process are shown. Moreover, the probability of survival and first passage time gives an estimation "
        "of the probability of failure depending on a chosen critical level $|u|$ that indicating malfunction."
        " It is important to note that the user-defined |u| in this case, defines the multiples of the standard deviation $\\sigma$. "
        "This means that $ |u|*\\sigma$ is passed to the functions. Further, the probability of survival is implemented by:"
        "$$L_{|X|}(u, t) = L_{|X|}(u, 0) \\exp \\left[-v_{|X|}^{+}(u) t\\right]\\:,$$"
        "where $ L_{|X|}(u, 0) $ is calculated with the CDF for the value $|u|$, and $ v_{|X|}^{+}(u) t $ describes the occurrence "
        "of crossing of level |u|."
        " The probability of failure time is given by:"
        "$$p_{T_{X}}(t) = v_{X}^{+}(u) \\exp \\left[-v_{|X|}^{+}(u) t\\right]\\:.$$"
        "For both functions, the Poisson approximation is used.",
        style={"font-size": "12px"},
        render_as_text=False,
        width=1200,
    )

    spacer = Div(text="", style={"font-size": "12px"}, width=300)

  
    return {
        "visu_or_theory_checkbox": visu_or_theory_checkbox,
        "covariance_checkbox": covariance_checkbox,
        "calculate_button_gauss": calculate_button_Gauss,
        "add_gaussian_button": Add_Gaussian_Button,
        "remove_gaussian_button": Remove_Gaussian_Button,
        "damping_slider": damping_slider,
        "stiffness_slider": stiffness_slider,
        "mass_slider": mass_slider,
        "sigma_slider": sigma_slider,
        "mean_slider": mean_slider,
        "correlation_slider": correlation_slider,
        "bandwidth_slider": bandwidth_slider,
        "frequency_slider": frequency_slider,
        "failure_slider": failure_slider,
        "title_div": title_div,
        "excitation_div": excitation_div,
        "parameter_div": para_div,
        "response_div": Response_div,
        "gaussian_div": Gaussian_div,
        "control_parameters_div": Control_para_div,
        "failure_div": Failure_div,
        "system_parameters_div": system_para_div,
        "sdof_system_div": SDOF_system_div,
        "excitation_parameters_div": excitation_para,
        "realization_div": realization_div,
        "frequency_title": frequency_title,
        "theory_div": theory_div,
        "spacer": spacer,
    }
    
    
def create_plots(standard_deviation):
 
    # excitation realization in time domain
    excitation_time = figure(
        title="Realization - Time Domain",
        x_axis_label="t [s]",
        y_axis_label="f(t)",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[0, 4],
        output_backend="svg",
    )

    style_figure(excitation_time)

    line_excitation_t = excitation_time.line([], [], line_color="royalblue", line_width=line_width)
    ds_excitation_t = line_excitation_t.data_source

    line_excitation_t_mean = excitation_time.line([], [], line_color="red", line_width=line_width)
    ds_excitation_t_mean = line_excitation_t_mean.data_source

    ex_new1 = excitation_time.line([], [], line_color="sienna", line_width=line_width)
    ex_new2 = excitation_time.line([], [], line_color="purple", line_width=line_width)
    ex_new3 = excitation_time.line([], [], line_color="green", line_width=line_width)
    ds_ex_new1 = ex_new1.data_source
    ds_ex_new2 = ex_new2.data_source
    ds_ex_new3 = ex_new3.data_source

    varea_excitation_time = excitation_time.varea(
        x=[], y1=[], y2=[], alpha=0.25, fill_color="olivedrab", legend_label="\u03C3"
    )
    ds_excitation_time_varea = varea_excitation_time.data_source

    varea_excitation_3_time = excitation_time.varea(
        x=[], y1=[], y2=[], alpha=0.25, fill_color="darkorange", legend_label="3\u03C3"
    )
    ds_excitation_time_3_varea = varea_excitation_3_time.data_source
    excitation_time.legend.label_text_font_size = "7pt"

    # excitation realization in frequency domain
    excitation_freq = figure(
        title="Realization - Frequency Domain",
        x_axis_label="f [Hz]",
        y_axis_label=" |F(f)| [Power/Hz]",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[0, 25],
        output_backend="svg",
    )
    style_figure(excitation_freq)

    line_excitation_f_real = excitation_freq.line([], [], line_color="royalblue", line_width=line_width)
    ds_excitation_f_real = line_excitation_f_real.data_source

    # covariance of the excitation
    covariance_exc = figure(
        title="Covariance",
        x_axis_label="\u03C4 [s]",
        y_axis_label="G_FF(\u03C4)",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[-4, 4],
        output_backend="svg",
    )
    style_figure(covariance_exc)

    line_covariance_exc = covariance_exc.line([], [], line_color="royalblue", line_width=line_width)
    ds_covariance_exc = line_covariance_exc.data_source

    spectral_density_exc = figure(
        title="Spectral Density",
        x_axis_label="f [Hz]",
        y_axis_label="|S_FF(f)|[Power/Hz]",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[0, 25],
        output_backend="svg",
    )
    style_figure(spectral_density_exc)
    spectral_density_exc.add_layout(
        LatexLabel(
            x=0,
            y=0,
            text="|S_XX(f)|",
            text_color="black",
            text_font_size="9pt",
            level="overlay",
            text_baseline="middle",
        )
    )

    line_spectral_density_exc = spectral_density_exc.line([], [], line_color="royalblue", line_width=line_width)
    ds_spectral_density_exc = line_spectral_density_exc.data_source

    # impulse response and transfer function
    transfer_time = figure(
        title="Impulse Response Function",
        y_axis_label="h_x(t)",
        x_axis_label="t [s]",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[-2, 2],
        output_backend="svg",
    )
    style_figure(transfer_time, scientific_y=True)

    line_transfer_t = transfer_time.line([], [], line_color="royalblue", line_width=line_width)
    ds_transfer_t = line_transfer_t.data_source

    transfer_freq_mag = figure(
        title="Harmonic Transfer Function",
        x_axis_label="f [Hz]",
        y_axis_label="|H_X(f)| [Power/Hz]",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[0, 25],
        output_backend="svg",
    )
    transfer_freq_mag.yaxis[0].formatter = PrintfTickFormatter(format="%4.0e")
    transfer_freq_mag.title.align = "center"
    transfer_freq_mag.toolbar.logo = None
    transfer_freq_mag.outline_line_width = line_width
    transfer_freq_mag.outline_line_color = "Black"
    transfer_freq_mag.title.text_font_size = text_font_size

    line_transfer_f_mag = transfer_freq_mag.line([], [], line_color="royalblue", line_width=line_width)
    ds_transfer_f_mag = line_transfer_f_mag.data_source

    # response realization in time domain
    response_time = figure(
        title="Realization - Time Domain",
        x_axis_label="t [s]",
        x_range=[0, 4],
        y_axis_label="x(t)",
        plot_width=plot_width,
        plot_height=plot_height,
        output_backend="svg",
    )
    style_figure(response_time, scientific_y=True)

    line_response_time = response_time.line([], [], line_color="royalblue", line_width=line_width)
    ds_response_time = line_response_time.data_source

    line_response_t_mean = response_time.line([], [], line_color="red", line_width=line_width)
    ds_response_t_mean = line_response_t_mean.data_source

    res_new = response_time.line([], [], line_color="sienna", line_width=line_width)
    res_new2 = response_time.line([], [], line_color="purple", line_width=line_width)
    res_new3 = response_time.line([], [], line_color="green", line_width=line_width)
    res_new_1 = res_new.data_source
    res_new_2 = res_new2.data_source
    res_new_3 = res_new3.data_source

    varea1_response = response_time.varea(
        x=[], y1=[], y2=[], alpha=0.25, fill_color="olivedrab", legend_label="\u03C3"
    )
    ds_response_varea1 = varea1_response.data_source

    varea2_response = response_time.varea(
        x=[], y1=[], y2=[], alpha=0.25, fill_color="darkorange", legend_label="3\u03C3"
    )
    ds_response_varea2 = varea2_response.data_source
    response_time.legend.label_text_font_size = "7pt"

    # response realization in frequency domain
    response_freq = figure(
        title="Realization - Frequency Domain",
        x_axis_label="f [Hz]",
        y_axis_label="|X(f)| [Power/Hz]",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[0, 25],
        output_backend="svg",
    )
    style_figure(response_freq, scientific_y=True)
    response_freq.title.text_font_size = text_font_size

    line_response_freq = response_freq.line([], [], line_color="royalblue", line_width=line_width)
    ds_responsef_freq = line_response_freq.data_source

    # covariance of the response
    covariance_res = figure(
        title="Covariance",
        x_axis_label="\u03C4 [s]",
        y_axis_label="G_XX(\u03C4)",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[-4, 4],
        output_backend="svg",
    )
    style_figure(covariance_res, scientific_y=True)

    line_covariance_res = covariance_res.line([], [], line_color="royalblue", line_width=line_width)
    ds_covariance_res = line_covariance_res.data_source

    line_covariance_res_proof = covariance_res.line([], [], line_color="green", line_width=line_width)
    ds_covariance_res_proof = line_covariance_res_proof.data_source

    # spectral density of the response
    spectral_density_res = figure(
        title="Spectral Density",
        x_axis_label="f [Hz]",
        y_axis_label="|S_XX(f)| [Power/Hz]",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[0, 25],
        output_backend="svg",
    )
    style_figure(spectral_density_res, scientific_y=True)
    spectral_density_res.title.text_font_size = text_font_size

    line_spectral_density_res = spectral_density_res.line([], [], line_color="royalblue", line_width=line_width)
    ds_spectral_density_res = line_spectral_density_res.data_source

    # probability density function of response
    pdf_gaussian = figure(
        title="Probability Density",
        x_axis_label="u",
        y_axis_label="p_X(u)",
        plot_width=plot_width,
        plot_height=plot_height,
        output_backend="svg",
    )
    style_figure(pdf_gaussian, scientific_y=True)

    line_gaussian_distribution = pdf_gaussian.line([], [], line_color="royalblue", line_width=line_width)
    ds_pdf_gaussian = line_gaussian_distribution.data_source

    vbar1_pdf_gaussian = pdf_gaussian.vbar(
        x=[],
        top=[],
        bottom=[],
        width=standard_deviation / 20,
        line_color="olivedrab",
        alpha=0.5,
        legend_label="\u03C3",
    )
    ds_pdf_gaussian_vbar1 = vbar1_pdf_gaussian.data_source

    vbar2_pdf_gaussian = pdf_gaussian.vbar(
        x=[],
        top=[],
        bottom=[],
        width=standard_deviation / 20,
        line_color="darkorange",
        alpha=0.5,
        legend_label="3\u03C3",
    )
    ds_pdf_gaussian_vbar2 = vbar2_pdf_gaussian.data_source
    pdf_gaussian.legend.label_text_font_size = "7pt"

    # cumulative distribution function of response
    cdf_gaussian = figure(
        title="Cumulative Distribution",
        x_axis_label="u",
        y_axis_label="\u03C6((u-\u03BC)/\u03C3)",
        plot_width=plot_width,
        plot_height=plot_height,
        output_backend="svg",
    )
    style_figure(cdf_gaussian, scientific_y=True)

    line_cdf_gaussian = cdf_gaussian.line([], [], line_color="royalblue", line_width=line_width)
    ds_cdf_gaussian = line_cdf_gaussian.data_source

    # probability of survival
    prob_of_survival = figure(
        title="Probability of Survival",
        x_axis_label="t [s]",
        y_axis_label="L_X(u = (u\u03C3),t)",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[0, 4],
        output_backend="svg",
    )
    style_figure(prob_of_survival)


    line_prob_of_survival = prob_of_survival.line([], [], line_color="royalblue", line_width=line_width)
    ds_prob_of_survival = line_prob_of_survival.data_source

    # first passage time
    prob_of_time = figure(
        title="First Passage Time",
        x_axis_label="t [s]",
        y_axis_label="p(T_X(t))",
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=[0, 4],
        output_backend="svg",
    )
    style_figure(prob_of_time, scientific_y=True)

    line_prob_of_time = prob_of_time.line([], [], line_color="royalblue", line_width=line_width)
    ds_prob_of_time = line_prob_of_time.data_source


    return {
        "figures": {
            "excitation_time": excitation_time,
            "excitation_freq": excitation_freq,
            "covariance_exc": covariance_exc,
            "spectral_density_exc": spectral_density_exc,
            "transfer_time": transfer_time,
            "transfer_freq_mag": transfer_freq_mag,
            "response_time": response_time,
            "response_freq": response_freq,
            "covariance_res": covariance_res,
            "spectral_density_res": spectral_density_res,
            "pdf_gaussian": pdf_gaussian,
            "cdf_gaussian": cdf_gaussian,
            "prob_of_survival": prob_of_survival,
            "prob_of_time": prob_of_time,
        },
        "sources": {
            "ds_excitation_t": ds_excitation_t,
            "ds_excitation_t_mean": ds_excitation_t_mean,
            "ds_ex_new1": ds_ex_new1,
            "ds_ex_new2": ds_ex_new2,
            "ds_ex_new3": ds_ex_new3,
            "ds_excitation_time_varea": ds_excitation_time_varea,
            "ds_excitation_time_3_varea": ds_excitation_time_3_varea,
            "ds_excitation_f_real": ds_excitation_f_real,
            "ds_covariance_exc": ds_covariance_exc,
            "ds_spectral_density_exc": ds_spectral_density_exc,
            "ds_transfer_t": ds_transfer_t,
            "ds_transfer_f_mag": ds_transfer_f_mag,
            "ds_response_time": ds_response_time,
            "ds_response_t_mean": ds_response_t_mean,
            "res_new_1": res_new_1,
            "res_new_2": res_new_2,
            "res_new_3": res_new_3,
            "ds_response_varea1": ds_response_varea1,
            "ds_response_varea2": ds_response_varea2,
            "ds_responsef_freq": ds_responsef_freq,
            "ds_covariance_res": ds_covariance_res,
            "ds_covariance_res_proof": ds_covariance_res_proof,
            "ds_spectral_density_res": ds_spectral_density_res,
            "ds_pdf_gaussian": ds_pdf_gaussian,
            "ds_pdf_gaussian_vbar1": ds_pdf_gaussian_vbar1,
            "ds_pdf_gaussian_vbar2": ds_pdf_gaussian_vbar2,
            "ds_cdf_gaussian": ds_cdf_gaussian,
            "ds_prob_of_survival": ds_prob_of_survival,
            "ds_prob_of_time": ds_prob_of_time,
        },
    }

def layout(widgets, figures):
    
    system_sliders = row(
    widgets["damping_slider"],
    widgets["stiffness_slider"],
    widgets["mass_slider"],
    widgets["frequency_title"],
    )

    transformation_plots = row(
    widgets["spacer"],
    figures["transfer_time"],
    figures["transfer_freq_mag"],
    widgets["spacer"],
    )

    covariance_sliders = row(
    widgets["sigma_slider"],
    widgets["mean_slider"],
    widgets["correlation_slider"],
    widgets["frequency_slider"],
    widgets["bandwidth_slider"],
    )

    realization_edit = row(
    widgets["spacer"],
    widgets["add_gaussian_button"],
    widgets["remove_gaussian_button"],
    widgets["spacer"],
    )

    excitation_plots = row(
    figures["covariance_exc"],
    figures["spectral_density_exc"],
    figures["excitation_time"],
    figures["excitation_freq"],
    )

    response_plots = row(
    figures["covariance_res"],
    figures["spectral_density_res"],
    figures["response_time"],
    figures["response_freq"],
    )

    gaussian_plots = row(
    widgets["spacer"],
    figures["pdf_gaussian"],
    figures["cdf_gaussian"],
    widgets["spacer"],
    )

    failure_plots = row(
    widgets["spacer"],
    figures["prob_of_survival"],
    figures["prob_of_time"],
    widgets["spacer"],
    )

    explain = widgets["theory_div"]

    visualization = column(
    widgets["control_parameters_div"],
    widgets["system_parameters_div"],
    system_sliders,
    widgets["excitation_parameters_div"],
    covariance_sliders,
    widgets["realization_div"],
    realization_edit,
    widgets["sdof_system_div"],
    transformation_plots,
    widgets["excitation_div"],
    excitation_plots,
    widgets["calculate_button_gauss"],
    widgets["response_div"],
    response_plots,
    widgets["gaussian_div"],
    gaussian_plots,
    widgets["failure_div"],
    row(widgets["spacer"], widgets["failure_slider"]),
    failure_plots,
    )

    final = column(
    widgets["title_div"],
    widgets["visu_or_theory_checkbox"],
    widgets["covariance_checkbox"],
    explain,
    visualization,
    )

    
    return explain, visualization, final