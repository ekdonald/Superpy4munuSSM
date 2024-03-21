Data files and limit/measurement implementation
===============================================

All limits in HiggsBounds and measurements in HiggsSignals are implemented
through `json` data files. The structure of these files is formally defined
through `json schemas <https://json-schema.org/>`_ that can be found
`here <https://gitlab.com/higgsbounds/higgstools/-/tree/develop/json>`_.

All HiggsBounds limit files should validate against
`LimitConditional.schema.json <https://gitlab.com/higgsbounds/higgstools/-/raw/develop/json/LimitConditional.schema.json>`_
(which is equivalent to
`Limit.schema.json <https://gitlab.com/higgsbounds/higgstools/-/raw/develop/json/Limit.schema.json>`_
but gives more useful error messages).
The following documentation is generated from
`Limit.schema.json <https://gitlab.com/higgsbounds/higgstools/-/raw/develop/json/Limit.schema.json>`_.


.. raw:: html
    :file: _static/Limit_schema.html
