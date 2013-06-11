%%
% AIAA class for technical conference papers
%
%  by Bil Kleb and Bill Wood
%
% TODO:
% - handout option?
% - submit option?

\NeedsTeXFormat{LaTeX2e}[1994/12/01]

\ProvidesClass{aiaa-tc}[2007/09/20 v4.0.1 AIAA technical conference paper]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}%
                \typeout{Unknown option: passing ``\CurrentOption"
                         option on to the underlying LaTeX article class.}}
\ProcessOptions

\LoadClass[letterpaper,twoside,10pt]{article}

%%
% Font

\RequirePackage{times}

%%
% Page layout

\RequirePackage[tmargin=0.25in,
                headheight=0.50in,
                headsep=0.25in,
                lmargin=1in,
                rmargin=1in,
                footskip=0.25in,
                bmargin=0.75in]{geometry}

%%
% Footnotes

\renewcommand\thefootnote{\@fnsymbol\c@footnote}

%%
% Page style (headers and footers)

\RequirePackage{fancyhdr}% simple control of headers and footers
\RequirePackage{lastpage}
\fancyhf{}
\renewcommand{\headrulewidth}{0pt}
\cfoot{\footnotesize\normalfont
       \thepage~of~\pageref{LastPage}\\
       \rule[.2\baselineskip]{0.5in}{0.2pt}\\
        American Institute of Aeronautics and Astronautics}
\pagestyle{fancy}

%%
% Title block -- original from article.cls

\renewcommand\maketitle{\par
  \begingroup
    \newpage
    \global\@topnum\z@
    \@maketitle
    \@thanks
  \endgroup
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{%
  \newpage
  \null
  \vskip 2em%
  \begin{center}%
  \let\footnote\thanks
    {\LARGE \@title \par}%
    \vskip 1.5em%
    {\large
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
  \end{center}%
  \par
  \vskip 1.5em}

%%
% Abstract environment

\renewenvironment{abstract}%
  {\small\quotation\bfseries}%
  {\normalfont\endquotation\normalsize}

%%
% Images

\RequirePackage{graphicx}
\setkeys{Gin}{ width=\linewidth, totalheight=\textheight, keepaspectratio }

%%
% Captions - original from latex.ltx

\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\footnotesize\bfseries #1\@. #2}%
  \ifdim \wd\@tempboxa >\hsize
    \footnotesize\bfseries #1\@. #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\setlength\belowcaptionskip{5pt plus 1pt minus 1pt}

%%
% Section headers - original from article.cls

\setcounter{secnumdepth}{3}

\renewcommand\thesection      {\@Roman\c@section.}
\renewcommand\thesubsection   {\@Alph\c@subsection.}
\renewcommand\thesubsubsection{\@arabic\c@subsubsection.}

\renewcommand\section{\@startsection {section}{1}{\z@}%
                                     {3.5ex \@plus -1ex \@minus -.2ex}%
                                     {2.3ex \@plus.2ex}%
                                     {\center\normalfont\large\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                       {-3.25ex\@plus -1ex \@minus -.2ex}%
                                       {1.5ex \@plus .2ex}%
                                       {\normalfont\normalsize\bfseries}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                       {-3.25ex\@plus -1ex \@minus -.2ex}%
                                       {1.5ex \@plus .2ex}%
                                       {\normalfont\normalsize\itshape}}
\renewcommand\paragraph{\@latex@error{Section level \noexpand\paragraph unavailable}\@eha}
\renewcommand\subparagraph{\@latex@error{Section level \noexpand\subparagraph unavailable}\@eha}

%%
% Colors

\RequirePackage{color}
  \definecolor{darkblue}{rgb}{0.0,0.0,0.5}

%%
% Hyperlinks

\RequirePackage[hyperfootnotes=false,
                hyperfigures,
                breaklinks,
                colorlinks,
                urlcolor=darkblue,
                filecolor=darkblue,
                citecolor=black,
                linkcolor=black,
                pagecolor=black,
                menucolor=black,
                anchorcolor=black]{hyperref}

%%
% References

\RequirePackage[numbers,sort&compress]{natbib}
\RequirePackage{hypernat}
  \bibpunct{}{}{,}{s}{}{} % superscript citation

% superscript reference list -- original from article.cls
\renewcommand{\@biblabel}[1]{\textsuperscript#1}
\renewenvironment{thebibliography}[1]%
     {\section*{\refname}
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\setlength{\leftmargin}{0pt}%
            \settowidth{\labelwidth}{\@biblabel{#1}}%
            \setlength{\itemindent}{\parindent}%
            \advance\itemindent by \labelwidth%
            \setlength{\labelsep}{0.0em}%
            \setlength{\itemsep}{-\smallskipamount}%
            \@openbib@code%
            \usecounter{enumiv}%
            \let\p@enumiv\@empty%
            \renewcommand\theenumiv{\@arabic\c@enumiv}%
            \footnotesize}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
      {\def\@noitemerr
        {\@latex@warning{Empty `thebibliography' environment}}%
     \endlist}

%%
% Floats -- make 'h' placement default

\def\fps@figure{htbp}
\def\fps@table{htbp}
