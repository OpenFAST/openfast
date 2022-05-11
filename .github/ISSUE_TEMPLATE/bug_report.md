---
name: Bug report
about: Report a bug to help us improve
title: 'Bug report'
labels: "Type: Bug"
---

<?
NOTE: Is your issue a bug or a usage question?
Usage questions are best answered in the Q&A category of
GitHub Discussions (https://github.com/openfast/openfast/discussions).

Also, consider how the NREL team will be able to understand your question with
no prior context. Help us help you by re-reading your
post and asking yourself whether someone else can reasonably
understand the question.

The following form is a template that should be completed entirely.
It uses GitHub's Markdown syntax (see https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).
When including text output from OpenFAST or lines of code, they should be
wrapped with three ticks (the other symbol on the ~ key in US keyboards)
before and after, like this:

```
**************************************************************************************************
 OpenFAST

 Copyright (C)  National Renewable Energy Laboratory
 Copyright (C)  Envision Energy USA LTD

 This program is licensed under Apache License Version 2.0 and comes with ABSOLUTELY NO WARRANTY.
 See the "LICENSE" file distributed with this software for details.
 **************************************************************************************************

 OpenFAST-v2.0.0
 Compile Info:
```
?>

# Acknowledgment of Due Diligence

I declare that I have done the following due diligence prior to creating this post
- [ ] Searched [GitHub Issues](https://github.com/openfast/openfast/issues) and [GitHub Discussions](https://github.com/openfast/openfast/discussions)
- [ ] Searched the [NREL Forum](https://forums.nrel.gov)
- [ ] Searched the internet (Google, Bing, etc)

# Bug description

<A clear and concise description of the bug.>

## To Reproduce

<Update the following list with your specific information.>

Steps to reproduce the behavior:
1. Compile with '...'
2. Run '...' case with '...' settings
3. Open '...' output
4. See the error

## Expected behavior

<A clear and concise description of what you expected to happen.>

## Text grabs or screenshots

<?
Add relevant text output or screenshots to help explain your problem.
Text output is better since it is searchable, but sometimes a screen
shot is more clear. Use your judgement. Please do not post a screenshot of text.
?>

## OpenFAST Version

<?
Please provide as much detail as possible including git commit.
The best information is the OpenFAST system description that prints when running OpenFAST:

```
**************************************************************************************************
 OpenFAST

 Copyright (C)  National Renewable Energy Laboratory
 Copyright (C)  Envision Energy USA LTD

 This program is licensed under Apache License Version 2.0 and comes with ABSOLUTELY NO WARRANTY.
 See the "LICENSE" file distributed with this software for details.
 **************************************************************************************************

 OpenFAST-v0.0.0
 Compile Info:
  - Architecture: 64 bit
  - Precision: double
  - Date: Nov 27 2018
  - Time: 17:19:38
 Execution Info:
  - Date: 11/29/2018
  - Time: 10:52:28-0700
```
?>

## System Information - Add your specific information:
 - OS: <e.g. Ubuntu 14.04 or macOS 10.12>
 - Compiler: <e.g. GFortran 4.4>
 - Compiler settings: <e.g. CMake flags or other settings>
