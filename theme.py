"""
theme.py - pastel turquoise styling and hand-built decorative SVG for app.py.

Everything here is drawn from scratch with <path> and <circle> elements, so
there are no image files, no external requests and no third-party licence to
honour. Change the palette in PALETTE and the whole app follows.
"""

import streamlit as st

#  Semantic names, so recolouring the app means editing only this dict.
#  'accent' drives the decorative squiggles and blobs; 'primary' drives the
#  buttons and the GC line on the chart.
PALETTE = {
    'bg':         '#FDF2F6',   # pale blush page
    'panel':      '#FBE3EC',   # sidebar / cards
    'ink':        '#5A2740',   # deep plum text
    'accent':     '#F2B8CE',   # mid pink, second helix strand
    'accent_dk':  '#E58FB4',   # stronger pink
    'primary':    '#35BE92',   # buttons, first helix strand
    'primary_dk': '#259B76',   # button hover
    'soft':       '#F7D3E0',   # hairlines, target band
    'btn_rest':   '#EAF7F2',   # button background at rest
}


# --------------------------------------------------------------------------
# Decorative SVG
# --------------------------------------------------------------------------
def helix(width: int = 1200, height: int = 70, period: int = 50,
          opacity: float = 0.75) -> str:
    """A DNA double helix drawn as two antiphase waves joined by rungs.

    Built from quadratic segments: one Q sets the first bump, then every T
    reflects the previous control point, so the wave alternates up/down on
    its own. The second strand uses a mirrored first control point, which
    puts it exactly half a period out of step.
    """
    mid = height / 2
    amp = height * 0.34
    n = width // period

    def strand(up_first: bool) -> str:
        c = mid - amp if up_first else mid + amp
        d = [f'M0,{mid}', f'Q{period/2},{c} {period},{mid}']
        d += [f'T{period * (i + 1)},{mid}' for i in range(1, n)]
        return ' '.join(d)

    # rungs sit at the quarter points, where the strands are furthest apart
    y_hi, y_lo = mid - amp * 0.75, mid + amp * 0.75
    rungs = ''.join(
        f'<line x1="{period/2 + i*period}" y1="{y_hi}" '
        f'x2="{period/2 + i*period}" y2="{y_lo}" '
        f'stroke="{PALETTE["accent"]}" stroke-width="2" stroke-linecap="round" '
        f'opacity="{0.35 + 0.3 * (i % 2)}"/>'
        for i in range(n)
    )

    return (
        f'<svg viewBox="0 0 {width} {height}" width="100%" height="{height}" '
        f'preserveAspectRatio="none" xmlns="http://www.w3.org/2000/svg" '
        f'role="presentation" style="opacity:{opacity}">'
        f'{rungs}'
        f'<path d="{strand(True)}" fill="none" stroke="{PALETTE["primary"]}" '
        f'stroke-width="3" stroke-linecap="round"/>'
        f'<path d="{strand(False)}" fill="none" stroke="{PALETTE["accent"]}" '
        f'stroke-width="3" stroke-linecap="round"/>'
        f'</svg>'
    )


def squiggle_rule(width: int = 1200, height: int = 24, period: int = 150) -> str:
    """A single loose wave, for use as a section divider.

    The period is deliberately long: the SVG is stretched to the container
    width with preserveAspectRatio="none", so a short period reads as a
    cramped zigzag once it is squeezed into a column.
    """
    mid = height / 2
    half = period / 2
    d = [f'M0,{mid}', f'Q{half},{mid - 8} {period},{mid}']
    d += [f'T{period * (i + 1)},{mid}' for i in range(1, width // period)]
    return (
        f'<svg viewBox="0 0 {width} {height}" width="100%" height="{height}" '
        f'preserveAspectRatio="none" xmlns="http://www.w3.org/2000/svg" '
        f'role="presentation">'
        f'<path d="{" ".join(d)}" fill="none" stroke="{PALETTE["soft"]}" '
        f'stroke-width="2.5" stroke-linecap="round"/></svg>'
    )


# --------------------------------------------------------------------------
# CSS
# --------------------------------------------------------------------------
def _css() -> str:
    p = PALETTE
    return f"""
    <style>
    /* floating "cells": filled blobs and membrane rings, fixed so they do
       not slide around while the results scroll */
    .stApp {{
        background-color: {p['bg']};
        background-image:
            radial-gradient(circle at 9% 16%,  rgba(242,184,206,.42) 0 54px, transparent 55px),
            radial-gradient(circle at 92% 9%,  rgba(229,143,180,.20) 0 86px, transparent 87px),
            radial-gradient(circle at 83% 64%, transparent 0 40px, rgba(242,184,206,.45) 41px 45px, transparent 46px),
            radial-gradient(circle at 16% 78%, transparent 0 26px, rgba(53,190,146,.22) 27px 30px, transparent 31px),
            radial-gradient(circle at 48% 94%, rgba(242,184,206,.30) 0 38px, transparent 39px);
        background-attachment: fixed;
        background-repeat: no-repeat;
    }}

    /* buttons: pale sage at rest, solid green on hover */
    .stButton > button, .stDownloadButton > button {{
        background: {p['btn_rest']};
        color: {p['ink']};
        border: 1.6px solid {p['accent']};
        border-radius: 12px;
        font-weight: 600;
        transition: all .15s ease;
    }}
    .stButton > button:hover, .stDownloadButton > button:hover {{
        background: {p['primary']};
        color: #FFFFFF;
        border-color: {p['primary_dk']};
        transform: translateY(-1px);
    }}
    .stButton > button[kind="primary"] {{
        background: {p['primary']};
        color: #FFFFFF;
        border-color: {p['primary_dk']};
    }}
    .stButton > button[kind="primary"]:hover {{
        background: {p['primary_dk']};
    }}

    /* rounded, tinted panels */
    [data-testid="stMetric"] {{
        background: rgba(255,255,255,.65);
        border: 1.5px solid {p['soft']};
        border-radius: 14px;
        padding: 12px 16px;
    }}
    .stTabs [data-baseweb="tab-list"] {{ gap: 6px; }}
    .stTabs [data-baseweb="tab"] {{
        background: rgba(255,255,255,.55);
        border-radius: 10px 10px 0 0;
        padding: 8px 16px;
    }}
    .stTabs [aria-selected="true"] {{
        background: {p['panel']};
        font-weight: 600;
    }}
    .stCode, pre {{ border-radius: 12px !important; }}
    textarea, .stTextArea textarea {{
        border-radius: 12px !important;
        border-color: {p['soft']} !important;
    }}
    h1 {{ letter-spacing: -.5px; }}
    </style>
    """


def apply_theme() -> None:
    """Inject the stylesheet. Call once, right after set_page_config."""
    st.markdown(_css(), unsafe_allow_html=True)


def draw(svg: str, height: int | None = None) -> None:
    """Render an SVG string.

    The SVG is base64-encoded into an <img> data URI rather than injected as
    raw markup: Streamlit's markdown sanitiser strips inline <svg> on some
    versions, and an <img> always survives. It also scopes any gradient or
    filter ids to the image, so two SVGs on one page cannot collide.
    """
    import base64
    b64 = base64.b64encode(svg.encode('utf-8')).decode('ascii')
    style = 'width:100%;display:block;'
    if height:
        style += f'height:{height}px;'
    st.markdown(
        f'<img src="data:image/svg+xml;base64,{b64}" style="{style}" alt=""/>',
        unsafe_allow_html=True,
    )
