\ProvidesClass{article-Bram}[2021/03/26 v1.2 article-Bram and Jelle]

% new package options
\newif\if@print
\DeclareOption{print}{\@printtrue}
\newif\if@serif
\DeclareOption{serif}{\@seriftrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

\ProcessOptions\relax
\LoadClass[a4paper, onecolumn, twoside]{article}

% Overall document layout
\RequirePackage[top=1in, bottom=1in, left=0.6in, right=0.6in, twoside]{geometry} % page setup
\RequirePackage[small]{titlesec} % titles layout
\titlespacing{\section}{0pt}{2ex}{1ex}
\titlespacing{\subsection}{10pt}{1.5ex}{0.75ex}
\titlespacing{\subsubsection}{10pt}{0.75ex}{0.5ex}
\setcounter{tocdepth}{3}
\RequirePackage{graphicx}\graphicspath{{./Figures/}}
\RequirePackage{parskip} % no indents and better paragraph rules
\usepackage{scrextend}\deffootnote[1em]{1em}{1em}{\textsuperscript{\thefootnotemark}\,} % footnote settings
\RequirePackage{fancyhdr} % fancy page layout, with header rule
\fancyhead[RO,LE]{\thepage}
\fancyhead[RE]{\leftmark}
\fancyhead[LO]{\rightmark}
\pagestyle{fancy}
\fancyfoot[c]{\thepage}
\usepackage{cancel}
% Fonts
\RequirePackage{amsmath} % needs before math font selection
\RequirePackage{fontspec} % advanced font selector
\setmainfont[Scale = MatchLowercase]{DejaVu Serif} % open serif font
\setsansfont[Scale = MatchLowercase]{DejaVu Sans} % open sans font
\setmonofont[Scale = MatchLowercase]{DejaVu Sans Mono}
\if@serif
\usepackage{unicode-math}
\setmathfont[Scale = MatchLowercase]{texgyredejavu-math.otf} % math font equal in style and size to DejaVu Serif
\else
\renewcommand{\familydefault}{\sfdefault}
\RequirePackage{arevmath} % math font that is equal in style DejaVu Sans, but better a, nu, omega, l, x ect
\DeclareMathSizes{8}{6.32}{5}{5}
\DeclareMathSizes{9}{7.1}{6}{5}
\DeclareMathSizes{10}{7.9}{6}{5}
\DeclareMathSizes{10.95}{8.65}{7}{5}
\DeclareMathSizes{12}{9.4}{7.6}{6}
\SetSymbolFont{symbols}{bold}{OMS}{zavm}{m}{n} % surpress warning of bold math sybols
\fi
\RequirePackage{bm} % bold greek and other symbols in math mode
\RequirePackage{csquotes} % proper quotation marks
\RequirePackage{physics} % propper fonts for derivatives, operators ect.

% Packages for extra enviroments/options
\RequirePackage[usenames, svgnames, table, xcdraw]{xcolor} % color selection by name
\RequirePackage[justification=centering]{caption} % caption settings
\captionsetup{font = small, labelsep = colon, labelfont=bf, belowskip=-10pt}
\RequirePackage{stfloats} % full size figures on the correct page
\RequirePackage[shortlabels]{enumitem} % list/items with advanced layout options
\setlist{nosep} % list/items with reduced white space
\RequirePackage{booktabs} % more table settings
\RequirePackage{longtable}
\RequirePackage{multirow} % tables with multirow cells
\RequirePackage[list=true]{subcaption} % subcaption settings
\DeclareCaptionSubType[alph]{figure} % subfigure caption alphabetic
\captionsetup[subfigure]{labelformat=brace,labelsep=space} % subcaption format
\newenvironment{code}{\captionsetup{type=listing}}{}
\RequirePackage{todonotes}

% Internal referencing and external links
\RequirePackage[hidelinks]{hyperref} % hyperrefs without ugly borders
\urlstyle{same} % urls in standard font
\if@print
\hypersetup{colorlinks=false} % color options for all links
\else
\hypersetup{colorlinks=true, urlcolor=MediumBlue, filecolor=MediumBlue, linkcolor=Crimson, citecolor=Green} % color options for links
\fi
\RequirePackage{cleveref} % smarter referencing to internal objects

% Math vector and matrices upright and bold
\def\vec#1{\ensuremath{\mathchoice
		{\mbox{\boldmath$\displaystyle\mathbf{#1}$}}
		{\mbox{\boldmath$\textstyle\mathbf{#1}$}}
		{\mbox{\boldmath$\scriptstyle\mathbf{#1}$}}
		{\mbox{\boldmath$\scriptscriptstyle\mathbf{#1}$}}}}

% Math exponential as operator e^
\let\exp\relax\DeclareMathOperator{\exp}{e}

% New commands
\renewcommand*\title[1]{\def\@title{#1}\hypersetup{pdftitle=#1}}
\newcommand*\subtitle[1]{\def\@subtitle{#1}}
%\newcommand*\keywords[1]{\def\@keywords{#1}\hypersetup{pdfkeywords={#1}}}
\newcommand*\institute[1]{\def\@institute{#1}}

% setup one column abstract
\RequirePackage{environ}
\RequirePackage{abstract}

\RenewEnviron{abstract}{%
	%\twocolumn[%
	\maketitle%
	\vspace{-0.2cm}
	\rule{\textwidth}{0.4pt}
	\vspace{-2\baselineskip}
	\section*{Abstract}%
	\BODY%
	
%	\ifdefined\@keywords
%		\textbf{Keywords:} \@keywords\\%
	\vspace{-0.3cm}
	\rule{\textwidth}{0.4pt}\\
	\vspace{-0.4cm}
	%]\saythanks%
}

% References
\RequirePackage[style=ieee, citestyle=numeric-comp, backend=biber]{biblatex} % IEEEtran cite
\addbibresource{bibliography.bib} % name and location of bibliography libary
\setlength\bibitemsep{0.3\baselineskip} % extra white space

\if@print
% Normal IEEEtran style with doi printed
\else
% Title as link to DOI reference
\ExecuteBibliographyOptions{doi=false}
\newbibmacro{string+doi}[1]{%
	\iffieldundef{doi}{#1}{\href{http://dx.doi.org/\thefield{doi}}{#1}}}
\DeclareFieldFormat{title}{\usebibmacro{string+doi}{\mkbibemph{#1}}}
\DeclareFieldFormat[article]{title}{\usebibmacro{string+doi}{\mkbibquote{#1}}}
\fi

\RequirePackage{siunitx}
\sisetup{separate-uncertainty=true, allow-number-unit-breaks, mode=text}

\definecolor{LightGray}{gray}{0.90}
\usepackage[newfloat]{minted}
\captionsetup[listing]{position=below}
\setminted{bgcolor=LightGray, linenos=true, xleftmargin=17pt, breaklines, breakanywhere, fontsize=\footnotesize}



% % dual column again
% \ProvidesClass{article-Bram}[2021/03/26 v1.2 article-Bram and Jelle]

% % new package options
% \newif\if@print
% \DeclareOption{print}{\@printtrue}
% \newif\if@serif
% \DeclareOption{serif}{\@seriftrue}
% \DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

% \ProcessOptions\relax
% \LoadClass[a4paper, twocolumn, twoside]{article}

% % Overall document layout
% \RequirePackage[top=1in, bottom=1in, left=0.6in, right=0.6in, twoside]{geometry} % page setup
% \RequirePackage[small]{titlesec} % titles layout
%     \titlespacing{\section}{0pt}{2ex}{1ex}
%     \titlespacing{\subsection}{10pt}{1.5ex}{0.75ex}
%     \titlespacing{\subsubsection}{10pt}{0.75ex}{0.5ex}
% \setcounter{tocdepth}{3}
% \RequirePackage{graphicx}\graphicspath{{./Figures/}}
% \RequirePackage{parskip} % no indents and better paragraph rules
% \usepackage{scrextend}\deffootnote[1em]{1em}{1em}{\textsuperscript{\thefootnotemark}\,} % footnote settings
% \RequirePackage{fancyhdr} % fancy page layout, with header rule
% \fancyhead[RO,LE]{\thepage}
% \fancyhead[RE]{\leftmark}
% \fancyhead[LO]{\rightmark}
% \pagestyle{fancy}
% \fancyfoot[c]{\thepage}
% \usepackage{cancel}
% % Fonts
% \RequirePackage{amsmath} % needs before math font selection
% \RequirePackage{fontspec} % advanced font selector
% \setmainfont[Scale = MatchLowercase]{DejaVu Serif} % open serif font
% \setsansfont[Scale = MatchLowercase]{DejaVu Sans} % open sans font
% \setmonofont[Scale = MatchLowercase]{DejaVu Sans Mono}
% \if@serif
% 	\usepackage{unicode-math}
% 	\setmathfont[Scale = MatchLowercase]{texgyredejavu-math.otf} % math font equal in style and size to DejaVu Serif
% \else
% \renewcommand{\familydefault}{\sfdefault}
%  	\RequirePackage{arevmath} % math font that is equal in style DejaVu Sans, but better a, nu, omega, l, x ect
% 	\DeclareMathSizes{8}{6.32}{5}{5}
% 	\DeclareMathSizes{9}{7.1}{6}{5}
% 	\DeclareMathSizes{10}{7.9}{6}{5}
% 	\DeclareMathSizes{10.95}{8.65}{7}{5}
% 	\DeclareMathSizes{12}{9.4}{7.6}{6}
%  	\SetSymbolFont{symbols}{bold}{OMS}{zavm}{m}{n} % surpress warning of bold math sybols
% \fi
% \RequirePackage{bm} % bold greek and other symbols in math mode
% \RequirePackage{minted}
% \RequirePackage{csquotes} % proper quotation marks
% \RequirePackage{physics} % propper fonts for derivatives, operators ect.

% % Packages for extra enviroments/options
% \RequirePackage[usenames, svgnames, table, xcdraw]{xcolor} % color selection by name
% \RequirePackage[justification=centering]{caption} % caption settings
% \captionsetup{font = small, labelsep = colon, labelfont=bf, belowskip=-10pt}
% \RequirePackage{stfloats} % full size figures on the correct page
% \RequirePackage[shortlabels]{enumitem} % list/items with advanced layout options
% \setlist{nosep} % list/items with reduced white space
% \RequirePackage{booktabs} % more table settings
% \RequirePackage{multirow} % tables with multirow cells
% \RequirePackage[list=true]{subcaption} % subcaption settings
% \DeclareCaptionSubType[alph]{figure} % subfigure caption alphabetic
% \captionsetup[subfigure]{labelformat=brace,labelsep=space} % subcaption format
% \newenvironment{code}{\captionsetup{type=listing}}{}
% \RequirePackage{todonotes}

% % Internal referencing and external links
% \RequirePackage[hidelinks]{hyperref} % hyperrefs without ugly borders
% \urlstyle{same} % urls in standard font
% \if@print
% 	\hypersetup{colorlinks=false} % color options for all links
% \else
% 	\hypersetup{colorlinks=true, urlcolor=MediumBlue, filecolor=MediumBlue, linkcolor=Crimson, citecolor=Green} % color options for links
% \fi
% \RequirePackage{cleveref} % smarter referencing to internal objects

% % Math vector and matrices upright and bold
% \def\vec#1{\ensuremath{\mathchoice
		% 		{\mbox{\boldmath$\displaystyle\mathbf{#1}$}}
		% 		{\mbox{\boldmath$\textstyle\mathbf{#1}$}}
		% 		{\mbox{\boldmath$\scriptstyle\mathbf{#1}$}}
		% 		{\mbox{\boldmath$\scriptscriptstyle\mathbf{#1}$}}}}

% % Math exponential as operator e^
% \let\exp\relax\DeclareMathOperator{\exp}{e}

% % New commands
% \renewcommand*\title[1]{\def\@title{#1}\hypersetup{pdftitle=#1}}
% \newcommand*\subtitle[1]{\def\@subtitle{#1}}
% \newcommand*\keywords[1]{\def\@keywords{#1}\hypersetup{pdfkeywords={#1}}}
% \newcommand*\institute[1]{\def\@institute{#1}}



% % setup one column abstract
% \RequirePackage{environ}
% \RequirePackage{abstract}

% \RenewEnviron{abstract}{%
	% 	\twocolumn[%
	% 		\maketitle%
	% 		\vspace{-0.2cm}
	% 		\rule{\textwidth}{0.4pt}
	% 		\vspace{-2\baselineskip}
	% 		\section*{Abstract}%
	% 			\BODY%
	
	% 		\ifdefined\@keywords
	% 			\textbf{Keywords:} \@keywords\\%
	% 			\vspace{-0.2cm}
	% 		\rule{\textwidth}{0.4pt}\\
	% 	]\saythanks%
	% }

% % References
% \RequirePackage[style=ieee, citestyle=numeric-comp, backend=biber]{biblatex} % IEEEtran cite
% \addbibresource{bibliography.bib} % name and location of bibliography libary
% \setlength\bibitemsep{0.3\baselineskip} % extra white space

% \if@print
% 	% Normal IEEEtran style with doi printed
% \else
% 	% Title as link to DOI reference
% 	\ExecuteBibliographyOptions{doi=false}
% 	\newbibmacro{string+doi}[1]{%
	% 	\iffieldundef{doi}{#1}{\href{http://dx.doi.org/\thefield{doi}}{#1}}}
% 	\DeclareFieldFormat{title}{\usebibmacro{string+doi}{\mkbibemph{#1}}}
% 	\DeclareFieldFormat[article]{title}{\usebibmacro{string+doi}{\mkbibquote{#1}}}
% \fi

% \RequirePackage{siunitx}
% \sisetup{separate-uncertainty=true, allow-number-unit-breaks}
