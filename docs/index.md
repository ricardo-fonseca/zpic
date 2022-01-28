---
title: zpic@edu

layout: splash

header:
  overlay_color: "#303030"
  overlay_filter: "0.05"
  overlay_image: /assets/images/twostream.png
  actions:
    - label: "Find us on GitHub"
      url: "https://github.com/ricardo-fonseca/zpic/"
#  caption: "PIC simulation of the Two-Stream instability"
excerpt: "Particle-in-Cell codes for plasma physics education and research"

intro: 
  - excerpt: 'The __ZPIC__ project is a suite of 1D/2D fully relativistic electromagnetic PIC codes, as well as 1D electrostatic. __ZPIC__ is geared towards plasma physics education and researchers looking for a simple, easily customizable, PIC code. Learn more about ZPIC [here](about).'
feature_row:
  - title: "Getting started"
    excerpt: "Learn the basics about ZPIC, including downloading and compiling instructions"
    url: "start"
    btn_label: "Start here"
    btn_class: "btn--primary"
  - title: "Documentation"
    excerpt: "Access user guides, algorithm details, and API reference"
    url: "documentation"
    btn_label: "Read the documentation"
    btn_class: "btn--primary"
  - title: "Examples"
    excerpt: "Check our library of well documented example ZPIC Jupyter notebooks"
    url: "examples"
    btn_label: "Notebook examples"
    btn_class: "btn--primary"
---

{% include feature_row id="intro" type="center"%}

{% include feature_row %}
