%-------------------------------------------------------------------------------
% Template class for ATLAS documents
%
% See ChangeLog or SVN log for history of changes.
% Responsible: Ian Brock (ian.brock@cern.ch)
%-------------------------------------------------------------------------------
\def\ATPackageVersion{01-03-02}

% --- Class structure: identification part
\ProvidesClass{latex/atlasdoc}[2014/09/22 v.01-03-02 ATLAS document class]
\NeedsTeXFormat{LaTeX2e}

% --- Class structure: declaration of options part
% Cover pages off and title page on by default
\newif\ifADOC@coverpage
\newif\ifADOC@preprint
\newif\ifADOC@auxmat
\newif\ifADOC@maketitle
\newif\ifADOC@aivpaper
\newif\ifADOC@letterpaper
\newif\ifADOC@CONFnote
\newif\ifADOC@PUBnote
\newif\ifADOC@koma
\newif\ifADOC@texmf
\ADOC@coverpagefalse
\ADOC@preprintfalse
\ADOC@auxmatfalse
\ADOC@maketitletrue
\ADOC@aivpapertrue
\ADOC@letterpaperfalse
\ADOC@CONFnotefalse
\ADOC@PUBnotefalse
\ADOC@komafalse
\ADOC@texmffalse
\DeclareOption{coverpage}{\ADOC@coverpagetrue}
% If preprint cover is made, regular title should be turned off
\DeclareOption{preprint}{\ADOC@preprinttrue \ADOC@auxmatfalse \ADOC@maketitlefalse}
% If auxiliary material cover is made, regular title should be turned off
\DeclareOption{auxmat}{\ADOC@preprintfalse \ADOC@auxmattrue \ADOC@maketitlefalse}
\DeclareOption{nomaketitle}{\ADOC@maketitlefalse}
% Paper size A4 by default
\DeclareOption{a4paper}{\ADOC@aivpapertrue \ADOC@letterpaperfalse}
\DeclareOption{letterpaper}{\ADOC@aivpaperfalse \ADOC@letterpapertrue}
% Type of document
\DeclareOption{CONF}{\ADOC@CONFnotetrue \ADOC@PUBnotefalse}
\DeclareOption{PUB}{\ADOC@CONFnotefalse \ADOC@PUBnotetrue}
% Use KOMA-script instead of standard LaTeX
\DeclareOption{koma}{\ADOC@komatrue}
% Packages available centrally, i.e. in texmf
\DeclareOption{texmf}{\ADOC@texmftrue}

% This class extends the article class
% Read all the atlasnote options, pass them to article
%\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

% --- Class structure: execution of options part
\ProcessOptions \relax

% --- Class structure: initial code part
% Defaults for title page
\newdimen\ADOC@skipbeforetitle
\ADOC@skipbeforetitle=60\p@
\newcommand{\skipbeforetitle}[1]{\ADOC@skipbeforetitle=#1}

\def\@title{Dummy title in atlasdoc\\ Use \textbackslash AtlasTitle to set it}
\newcommand{\AtlasTitleText}{Dummy AtlasTitleText in atlasdoc\\ Use \textbackslash AtlasTitle to set it}
\newcommand{\AtlasTitle}[1]{%
  \def\@title{#1}%
  \renewcommand{\AtlasTitleText}{#1}%
}
\def\@author{The ATLAS Collaboration}
\def\@date{\today}

% Author(s)
%\def\@makefnmark{\hbox{$^{\@thefnmark)}$}}
%\def\author#1{%% Treat the list of authors
%\setcounter{footnote}{0}\def\@currentlabel{}%
%\begingroup\def\thefootnote{\arabic{footnote}}
%\global\@topnum\z@ \large\begin{center}{\lineskip.5em
%\begin{tabular}[t]{c}#1\end{tabular}\par}
%\end{center}\par\vskip1.5em\@thanks\endgroup}
%\newenvironment{Authlist}{\begin{center}}{\end{center}}

% Document title
%\newcommand{\AtlasTitle}[1]%
%{
%  \def\@title{#1}
%  \def\@atlastitle{#1}
%}

% Draft heading code
\newcommand{\olddraftheading}{%
  % 
  % compute the time in hours and minutes
  % make new variables \timehh and \timemm
  % 
  \newcount\timehh\newcount\timemm
  \timehh=\time 
  \divide\timehh by 60 \timemm=\time
  \count255=\timehh\multiply\count255 by -60 \advance\timemm by \count255
  % put in header
  \markboth%
  % --- Left mark
  {\hspace*{\fill}{\protect\small\bf \fbox{DRAFT}}
    \hspace*{\fill}
    \protect\makebox[0pt][r]{\protect\small\sl\today\ -- 
      \ifnum\timehh<10 0\fi\number\timehh\,:\,\ifnum\timemm<10 0\fi\number\timemm}
  }%--- right mark
  {\protect\makebox[0pt][l]{\protect\small\sl\today\ -- 
      \ifnum\timehh<10 0\fi\number\timehh\,:\,\ifnum\timemm<10 0\fi\number\timemm}
    \hspace*{\fill}
    {\fbox{\textbf{DRAFT}}}
    \hspace*{\fill}}
  %
  \pagestyle{myheadings}
}

\newcommand{\draftheading}{%
  \RequirePackage{fancyhdr}
  %
  % Compute the time in hours and minutes; make new variables \timehh and \timemm
  %
  \newcount\timehh\newcount\timemm
  \timehh=\time 
  \divide\timehh by 60 \timemm=\time
  \count255=\timehh\multiply\count255 by -60 \advance\timemm by \count255
  %
  % Now make a string with date and time
  %
  \newcommand{\dateandtime}{\today \ --
    \ifnum\timehh<10 0\fi\number\timehh:\ifnum\timemm<10 0\fi\number\timemm}
  % \RequirePackage{scrdate,scrtime}
  \chead[]{\protect\small\fbox{DRAFT}}
  \rhead{}
  % \lhead[]{\ISOToday{} \thistime}
  % \lhead[]{\today{} \thistime}
  %\rhead{\thepage}
  \lfoot{\dateandtime}
  \cfoot{}
  \rfoot{\thepage}
  \setlength\headheight{16pt}
  \pagestyle{fancy}
}

% Draft version 
\def\ADOC@draftversion{}
\newcommand*{\draftversion}[1]{\def\ADOC@draftversion{#1}}
\newcommand*{\AtlasVersion}[1]{\def\ADOC@draftversion{#1}}

% ATLAS reference code - empty by default
% Also create variable with web page
\def\ADOC@refcode{}
\newcommand*{\AtlasRefCode}[1]{\def\ADOC@refcode{#1}}
\newcommand*{\AtlasPaperWebpage}{https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/\ADOC@refcode}

% ATLAS note number
\def\ADOC@atlasnote{}
\newcommand*{\AtlasNote}[1]{\def\ADOC@atlasnote{#1}}

% Journal
\def\ADOC@journal{}
\newcommand*{\journal}[1]{\def\ADOC@journal{#1}}
\newcommand*{\AtlasJournal}[1]{\def\ADOC@journal{#1}}

% Preprint
\def\ADOC@preprintnumber{}
\newcommand*{\PreprintIdNumber}[1]{\def\ADOC@preprintnumber{#1}}

% arXiv
\def\AtlasArXivId{}
\def\AtlasArXivURL{}
\newcommand*{\arXivId}[1]{\def\AtlasArXivId{#1}\def\AtlasArXivURL{http://arXiv.org/abs/#1}}
\newcommand*{\AtlasArXivCat}{[hep-ex]}

% HEP Data
\def\AtlasHepDataRecord{}
\newcommand*{\HepDataRecord}[1]{\def\AtlasHepDataRecord{http://hepdata.cedar.ac.uk/view/ins#1}}

% Submission details
\def\ADOC@subprepend{To be submitted to }
\newcommand*{\subprepend}[1]{\def\ADOC@subprepend{#1}}

% Abstract text
\def\AtlasAbstractText{}
\newcommand{\AtlasAbstract}[1]{\def\AtlasAbstractText{#1}}
\newcommand{\abstracttext}[1]{\def\AtlasAbstractText{#1}}

% Shortcut to add submission text
\newcommand{\submitted}{%
  \ifx\ADOC@journal\@empty
    \mbox{}%
  \else
    \vskip 0.8in plus 2in%
    \begin{center}%
      \ADOC@subprepend\space \textit{\ADOC@journal}%
      \vskip 1.5em%
    \end{center}%
    \vskip 0.5in
  \fi
}

% --- Class structure: invoke main class
\ifADOC@koma%
  \ifADOC@aivpaper%
    \LoadClass[titlepage,a4paper,11pt]{scrartcl}
  \else
    \LoadClass[titlepage,letterpaper,11pt]{scrartcl}
  \fi
\else
  \ifADOC@aivpaper%
    \LoadClass[titlepage,a4paper,11pt]{article}
  \else
    \LoadClass[titlepage,letterpaper,11pt]{article}
  \fi
\fi

\RequirePackage{geometry}

% Redefine the page size - use fraction of paper, so it works for A4 and letter
\geometry{
  centering,scale=0.75,%
  marginparsep=0.2cm,% 
  marginparwidth=1.8cm%
}

% Turn on bold math when using \bfseries
\def\bfseries{\fontseries\bfdefault\selectfont\boldmath}

% --- Class structure: main code part
% ATLAS draft cover page
\ifADOC@coverpage%
  \ifADOC@CONFnote
    \ifADOC@texmf
      \RequirePackage[CONF,atlasdoc]{atlascover}
    \else
      \RequirePackage[CONF,atlasdoc]{latex/atlascover}
    \fi
  \else
    \ifADOC@PUBnote
      \ifADOC@texmf
        \RequirePackage[PUB,atlasdoc]{atlascover}
      \else
        \RequirePackage[PUB,atlasdoc]{latex/atlascover}
      \fi
    \else
      \ifADOC@texmf
        \RequirePackage[atlasdoc]{atlascover}
      \else
        \RequirePackage[atlasdoc]{latex/atlascover}
      \fi
    \fi
  \fi
\fi
% CERN preprint cover page
\ifADOC@preprint%
  \ifADOC@texmf
    \RequirePackage[atlasdoc]{atlascover}
  \else
    \RequirePackage[atlasdoc]{latex/atlascover}
  \fi
\fi
% Auxiliary material cover page
\ifADOC@auxmat%
  \ifADOC@texmf
    \RequirePackage[atlasdoc]{atlascover}
  \else
    \RequirePackage[atlasdoc]{latex/atlascover}
  \fi
\fi

% Remove the numbers from all the headings (\section, \subsection, etc.)
%\setcounter{secnumdepth}{-1}

% Required packages
\RequirePackage{graphicx}
\RequirePackage{xcolor}
\RequirePackage{babel}
\RequirePackage[varg]{txfonts}
%\RequirePackage{mathptmx}
% Add ability to include line numbers
\RequirePackage{lineno}
% If you have an old lineno or do not have it, you can get it from the latex/other directory
% \RequirePackage{latex/other/lineno}


% Environment for analysis contributers
\newenvironment{atlascontribute}[1][3.5cm]
  {\section*{List of contributions to analysis}
   \begin{list}{}%
    {\setlength{\leftmargin}{#1}
     \setlength{\labelwidth}{#1}
     \addtolength{\leftmargin}{0.2cm}
     \setlength{\labelsep}{0.2cm}
     \setlength{\itemsep}{0.0cm}
     \renewcommand{\makelabel}[1]{##1 \hfill}}
  }
  {\end{list}
}

% Title page design
\newcommand{\maketitlepage}{%
  \null
  \thispagestyle{empty}
  \pagestyle{plain}
  \pagenumbering{arabic}
  \ifx\ADOC@draftversion\@empty\mbox{}%
  \else \flushright Draft version \ADOC@draftversion \\ 
  \fi
  \begin{center}
    \vspace*{-6mm}
    \includegraphics[width=3cm]{AT_atlaslogo} 
    \hfill 
    \begin{minipage}[b]{7cm}
      \vfill
      \begin{center}
        {\huge \bf \textsf{ATLAS NOTE}}%
        % ATLAS note number, if any 
        \ifx\ADOC@atlasnote\@empty
          % ATLAS reference code, if any and no note number given
          \ifx\ADOC@refcode\@empty \vskip 4mm \mbox{}%
          \else \vskip 4mm {\Large \ADOC@refcode} \vskip 0mm%
          \fi
        \else \vskip 4mm {\Large \ADOC@atlasnote} \vskip 0mm%
        \fi
        % Date, if any
        \ifx\@date\@empty\mbox{}%	
        \else \vskip 2mm \@date%
        \fi
      \end{center}
      \vfill
    \end{minipage}
    \hfill
    \includegraphics[width=3cm]{AT_cernlogo}
  \end{center}
  % Add line numbers to front page content 
  {%
    \ifx\ADOC@draftversion\@empty%
    \else\linenumbers
    \fi
    % Title, author(s)
    %\vskip 60\p@
    \vskip \ADOC@skipbeforetitle
    \begin{center}%
      {\Large\bf \@title \par}%
      \vskip 3em%
      {\large
       \lineskip .75em%
       \begin{tabular}[t]{c}%
         \@author
       \end{tabular}\par}%
    \end{center}\par
    \vfill
    % Abstract
    \if\AtlasAbstractText\@empty\mbox{}%
    \else%
      \begin{center}
        {\bfseries \abstractname}
        \quotation
        \setlength{\parindent}{0pt}\setlength{\parskip}{0.5ex}\noindent\AtlasAbstractText
        \endquotation
      \end{center}
    \fi
    \vfill
    % Submission details
    \submitted
    % Thanks
    \@thanks
  }%
}

% This some before-and-after code that surrounds the title page. It
% shouldn't need to be modified. I've pulled out the part the actually
% typesets the title page and placed it in the \maketitlepage command
% above.
\renewcommand\maketitle{%
  \begin{titlepage}%
    \let\footnotesize\small%
    \let\footnoterule\relax%
    \let \footnote \thanks%
    \maketitlepage%
    \let\thefootnote\relax\footnotetext{\copyright\ Copyright \the\year\ CERN for the benefit of the ATLAS Collaboration. \\ \indent \indent Reproduction of this article or parts of it is allowed as specified in the CC-BY-3.0 license.}
  \end{titlepage}%
  \setcounter{footnote}{0}%
  \setcounter{page}{1}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}


\AtBeginDocument{
  % If selected, generate ATLAS Cover Page
  \ifADOC@coverpage%
    \ACOV@PrintAtlasCover%
    \ACOV@PrintAtlasPage%
  \fi
  % If selected, generate ATLAS preprint front page
  \ifADOC@preprint%
    \ACOV@PrintPreprintCover%
  \fi
  % If selected, generate ATLAS auxiliary material front page
  \ifADOC@auxmat%
    \ACOV@PrintAuxmatCover%
  \fi
  % Automatically generate title page
  \ifADOC@maketitle
    \maketitle
  \fi
  % Add draft heading and line numbers
  \ifx\ADOC@draftversion\@empty%
  \else \linenumbers \draftheading 
  \fi
  % Hack for tikz (must be loaded here in order to avoid conflicts with graphicx user in \maketitle}
  % This hack does not appear to be necessary any more (ICB: 2014-08-26)
  %\if@usetikz \RequirePackage{tikz} \fi
  %\pagenumbering{arabic} % reset page numbering after title page
}

\endinput
