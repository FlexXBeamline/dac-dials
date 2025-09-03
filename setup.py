from setuptools import setup

setup(
    name="dac-dials",
    packages=["dac_dials"],
    version="0.1.0",
    description="DIALS add-on for MX in a diamond anvil cell (DAC)",
    author="Steve P. Meisburger",
    author_email="spm82@cornell.edu",
    url="https://github.com/FlexXBeamline/dac-dials",
    license="BSD",
    python_requires=">=3.10",
    install_requires=[
        "dials",
    ],
    entry_points={
        "console_scripts": [
            "dac.fit_gasket_mask=dac_dials.fit_gasket_mask_v2:run",
            "dac.fit_gasket_mask_v0=dac_dials.fit_gasket_mask:run",
            "dac.trim_image_range=dac_dials.trim_image_range:run",
        ],
        "libtbx.dispatcher.script": [
            "dac.fit_gasket_mask=dac_dials.fit_gasket_mask_v2:run",
            "dac.fit_gasket_mask_v0=dac_dials.fit_gasket_mask:run",
            "dac.trim_image_range=dac_dials.trim_image_range:run",
        ],
        "libtbx.precommit": ["dac_dials=dac_dials"],  # necessary?
    },
    include_package_data=True,
)
