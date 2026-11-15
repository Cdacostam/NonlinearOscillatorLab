# Verbatim copies of the LaTeX files of Supplementary Appendix B (revision v5, 7 Nov 2026).
# Do not edit: verify_sextic.py --tex DIR checks that these strings equal the .tex sources.
APPB_PHIPSI = r"""\begin{align}
\Phi(\lambda,\mu)={}&  b + 4a\lambda + 2x_0^{2}b\lambda - 12C_0\lambda^{2} + 4x_0^{2}a\lambda^{2} + x_0^{4}b\lambda^{2}\nonumber\\
&\quad  - 16C_0\mu + 2x_0^{4}b\mu + 4x_0^{4}a\lambda \mu + 2x_0^{6}b\lambda \mu + x_0^{8}b\mu^{2}\label{eq:Phi}
\end{align}
\begin{align}
\Psi(\lambda,\mu)={}&  2c + 3b\lambda + 6x_0^{2}c\lambda + 6a\lambda^{2} + 6x_0^{2}b\lambda^{2} + 6x_0^{4}c\lambda^{2}\nonumber\\
&\quad  - 12C_0\lambda^{3} + 6x_0^{2}a\lambda^{3} + 3x_0^{4}b\lambda^{3} + 2x_0^{6}c\lambda^{3} + 24a\mu + 6x_0^{4}c\mu\nonumber\\
&\quad  - 144C_0\lambda \mu + 24x_0^{2}a\lambda \mu + 6x_0^{4}b\lambda \mu + 12x_0^{6}c\lambda \mu + 6x_0^{4}a\lambda^{2} \mu + 6x_0^{6}b\lambda^{2} \mu\nonumber\\
&\quad  + 6x_0^{8}c\lambda^{2} \mu + 24x_0^{4}a\mu^{2} + 6x_0^{8}c\mu^{2} + 3x_0^{8}b\lambda \mu^{2} + 6x_0^{10}c\lambda \mu^{2} + 2x_0^{12}c\mu^{3}\label{eq:Psi}
\end{align}
"""

APPB_SEXTIC = r"""\begin{align}
e_0={}&  576\alpha^{2} \beta^{2} + 1\,536C_0\alpha \beta \gamma + 1\,024C_0^{2}\gamma^{2}\nonumber\\
e_1={}&  2\,304\alpha^{3} \beta + 1\,728x_0^{2}\alpha^{2} \beta^{2} - 2\,304C_0\alpha \beta^{2} + 3\,072C_0\alpha^{2} \gamma + 4\,608C_0x_0^{2}\alpha \beta \gamma\nonumber\\
&\quad  - 3\,072C_0^{2}\beta \gamma + 3\,072C_0^{2}x_0^{2}\gamma^{2}\nonumber\\
e_2={}&  4\,608x_0^{2}\alpha^{3} \beta - 18\,432C_0\alpha^{2} \beta + 2\,160x_0^{4}\alpha^{2} \beta^{2} - 4\,608C_0x_0^{2}\alpha \beta^{2} + 576C_0x_0^{4}\beta^{3}\nonumber\\
&\quad  - 768x_0^{4}\alpha^{3} \gamma + 6\,144C_0x_0^{2}\alpha^{2} \gamma - 24\,576C_0^{2}\alpha \gamma + 4\,992C_0x_0^{4}\alpha \beta \gamma - 6\,144C_0^{2}x_0^{2}\beta \gamma\nonumber\\
&\quad  + 768C_0^{2}x_0^{4}\gamma^{2}\nonumber\\
e_3={}&  3\,456x_0^{4}\alpha^{3} \beta - 18\,432C_0x_0^{2}\alpha^{2} \beta + 36\,864C_0^{2}\alpha \beta + 1\,440x_0^{6}\alpha^{2} \beta^{2} - 1\,152C_0x_0^{4}\alpha \beta^{2}\nonumber\\
&\quad  + 1\,152C_0x_0^{6}\beta^{3} - 1\,536x_0^{6}\alpha^{3} \gamma + 10\,752C_0x_0^{4}\alpha^{2} \gamma - 24\,576C_0^{2}x_0^{2}\alpha \gamma + 49\,152C_0^{3}\gamma\nonumber\\
&\quad  + 2\,304C_0x_0^{6}\alpha \beta \gamma - 7\,680C_0^{2}x_0^{4}\beta \gamma - 3\,584C_0^{2}x_0^{6}\gamma^{2}\nonumber\\
e_4={}&  1\,152x_0^{6}\alpha^{3} \beta - 4\,608C_0x_0^{4}\alpha^{2} \beta + 540x_0^{8}\alpha^{2} \beta^{2} + 1\,152C_0x_0^{6}\alpha \beta^{2} - 9\,216C_0^{2}x_0^{4}\beta^{2}\nonumber\\
&\quad  + 864C_0x_0^{8}\beta^{3} - 1\,152x_0^{8}\alpha^{3} \gamma + 7\,680C_0x_0^{6}\alpha^{2} \gamma - 18\,432C_0^{2}x_0^{4}\alpha \gamma - 480C_0x_0^{8}\alpha \beta \gamma\nonumber\\
&\quad  - 4\,608C_0^{2}x_0^{6}\beta \gamma - 576C_0^{2}x_0^{8}\gamma^{2}\nonumber\\
e_5={}&  144x_0^{8}\alpha^{3} \beta + 108x_0^{10}\alpha^{2} \beta^{2} + 432C_0x_0^{8}\alpha \beta^{2} + 288C_0x_0^{10}\beta^{3} - 384x_0^{10}\alpha^{3} \gamma\nonumber\\
&\quad  + 1\,728C_0x_0^{8}\alpha^{2} \gamma - 864C_0x_0^{10}\alpha \beta \gamma + 5\,184C_0^{2}x_0^{8}\beta \gamma + 1\,728C_0^{2}x_0^{10}\gamma^{2}\nonumber\\
e_6={}&  9x_0^{12}\alpha^{2} \beta^{2} + 36C_0x_0^{12}\beta^{3} - 48x_0^{12}\alpha^{3} \gamma - 216C_0x_0^{12}\alpha \beta \gamma - 432C_0^{2}x_0^{12}\gamma^{2}\label{eq:ecoef}
\end{align}
"""

APPB_G = r"""- 3\alpha^{2} \beta^{2} - 12C_0\beta^{3} + 16\alpha^{3} \gamma + 72C_0\alpha \beta \gamma + 144C_0^{2}\gamma^{2}
"""

APPB_W = r"""\begin{align}
W={}& - 324x_0^{4}\alpha^{6} \beta^{4} - 1\,296C_0x_0^{2}\alpha^{5} \beta^{4} + 3\,888C_0^{2}\alpha^{4} \beta^{4} - 2\,592C_0x_0^{4}\alpha^{4} \beta^{5}\nonumber\\
&\quad  - 7\,776C_0^{2}x_0^{2}\alpha^{3} \beta^{5} + 20\,736C_0^{3}\alpha^{2} \beta^{5} - 324C_0x_0^{6}\alpha^{3} \beta^{6} - 5\,832C_0^{2}x_0^{4}\alpha^{2} \beta^{6}\nonumber\\
&\quad  - 648C_0^{2}x_0^{6}\alpha \beta^{7} - 81C_0^{2}x_0^{8}\beta^{8} - 864x_0^{6}\alpha^{6} \beta^{3} \gamma - 1\,728C_0x_0^{4}\alpha^{5} \beta^{3} \gamma\nonumber\\
&\quad  + 17\,280C_0^{2}x_0^{2}\alpha^{4} \beta^{3} \gamma - 20\,736C_0^{3}\alpha^{3} \beta^{3} \gamma - 4\,320C_0x_0^{6}\alpha^{4} \beta^{4} \gamma - 10\,368C_0^{2}x_0^{4}\alpha^{3} \beta^{4} \gamma\nonumber\\
&\quad  + 10\,368C_0^{3}x_0^{2}\alpha^{2} \beta^{4} \gamma + 55\,296C_0^{4}\alpha \beta^{4} \gamma - 648C_0x_0^{8}\alpha^{3} \beta^{5} \gamma - 9\,072C_0^{2}x_0^{6}\alpha^{2} \beta^{5} \gamma\nonumber\\
&\quad  - 12\,096C_0^{3}x_0^{4}\alpha \beta^{5} \gamma - 864C_0^{2}x_0^{8}\alpha \beta^{6} \gamma - 864C_0^{3}x_0^{6}\beta^{6} \gamma - 108C_0^{2}x_0^{10}\beta^{7} \gamma\nonumber\\
&\quad  - 864x_0^{8}\alpha^{6} \beta^{2} \gamma^{2} - 2\,304C_0x_0^{6}\alpha^{5} \beta^{2} \gamma^{2} + 16\,128C_0^{2}x_0^{4}\alpha^{4} \beta^{2} \gamma^{2} - 9\,216C_0^{3}x_0^{2}\alpha^{3} \beta^{2} \gamma^{2}\nonumber\\
&\quad  - 13\,824C_0^{4}\alpha^{2} \beta^{2} \gamma^{2} - 2\,304C_0x_0^{8}\alpha^{4} \beta^{3} \gamma^{2} - 5\,760C_0^{2}x_0^{6}\alpha^{3} \beta^{3} \gamma^{2} - 23\,040C_0^{3}x_0^{4}\alpha^{2} \beta^{3} \gamma^{2}\nonumber\\
&\quad  + 69\,120C_0^{4}x_0^{2}\alpha \beta^{3} \gamma^{2} + 36\,864C_0^{5}\beta^{3} \gamma^{2} - 432C_0x_0^{10}\alpha^{3} \beta^{4} \gamma^{2} - 4\,320C_0^{2}x_0^{8}\alpha^{2} \beta^{4} \gamma^{2}\nonumber\\
&\quad  - 19\,008C_0^{3}x_0^{6}\alpha \beta^{4} \gamma^{2} - 5\,760C_0^{4}x_0^{4}\beta^{4} \gamma^{2} - 288C_0^{2}x_0^{10}\alpha \beta^{5} \gamma^{2} - 1\,152C_0^{3}x_0^{8}\beta^{5} \gamma^{2}\nonumber\\
&\quad  - 36C_0^{2}x_0^{12}\beta^{6} \gamma^{2} - 384x_0^{10}\alpha^{6} \beta \gamma^{3} - 2\,304C_0x_0^{8}\alpha^{5} \beta \gamma^{3} + 9\,216C_0^{2}x_0^{6}\alpha^{4} \beta \gamma^{3}\nonumber\\
&\quad  + 30\,720C_0^{3}x_0^{4}\alpha^{3} \beta \gamma^{3} - 129\,024C_0^{4}x_0^{2}\alpha^{2} \beta \gamma^{3} + 110\,592C_0^{5}\alpha \beta \gamma^{3} - 384C_0x_0^{10}\alpha^{4} \beta^{2} \gamma^{3}\nonumber\\
&\quad  - 1\,536C_0^{2}x_0^{8}\alpha^{3} \beta^{2} \gamma^{3} - 18\,432C_0^{4}x_0^{4}\alpha \beta^{2} \gamma^{3} + 55\,296C_0^{5}x_0^{2}\beta^{2} \gamma^{3} - 96C_0x_0^{12}\alpha^{3} \beta^{3} \gamma^{3}\nonumber\\
&\quad  - 576C_0^{2}x_0^{10}\alpha^{2} \beta^{3} \gamma^{3} - 5\,760C_0^{3}x_0^{8}\alpha \beta^{3} \gamma^{3} - 9\,984C_0^{4}x_0^{6}\beta^{3} \gamma^{3} - 64x_0^{12}\alpha^{6} \gamma^{4}\nonumber\\
&\quad  - 768C_0x_0^{10}\alpha^{5} \gamma^{4} - 768C_0^{2}x_0^{8}\alpha^{4} \gamma^{4} + 14\,336C_0^{3}x_0^{6}\alpha^{3} \gamma^{4} + 9\,216C_0^{4}x_0^{4}\alpha^{2} \gamma^{4}\nonumber\\
&\quad  - 110\,592C_0^{5}x_0^{2}\alpha \gamma^{4} + 110\,592C_0^{6}\gamma^{4}\nonumber
\end{align}
"""

