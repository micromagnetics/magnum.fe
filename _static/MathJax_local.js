document.write('<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML"></script>');
document.write('\
<script type="text/x-mathjax-config">\
          MathJax.Hub.Config({\
            "TeX": {\
              "Macros": {\
                "vec": ["{\\\\boldsymbol #1}", 1],\
                "mat": ["{\\\\boldsymbol #1}", 1],\
                "dx":  "{\\\\;\\\\text{d}\\\\vec{x}}"\
              }\
            },\
            "HTML-CSS": {\
              "linebreaks": {"automatic": true, "width": "container"}\
            }\
          });\
        </script>');
