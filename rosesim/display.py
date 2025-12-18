
import copy

import matplotlib
import matplotlib.pyplot as plt

from astropy import visualization as aviz

from mpl_toolkits.axes_grid1 import make_axes_locatable


IMG_CMAP = copy.copy(matplotlib.cm.get_cmap("viridis"))
IMG_CMAP.set_bad(color="black")


def show_image_wcs(
    image,
    percl=99,
    percu=None,
    vmin=None,
    vmax=None,
    is_mask=False,
    figsize=(10, 10),
    cmap="Greys",
    log=False,
    clip=True,
    show_colorbar=False,
    show_ticks=False,
    fig=None,
    ax=None,
    input_ratio=None,
    scale_args=None,
):
    """
    Show an image in matplotlib with some basic astronomically-appropriate stretching.
    From https://github.com/astropy/ccd-reduction-and-photometry-guide/blob/main/notebooks/convenience_functions.py

    Parameters
    ----------
    image (numpy array): The image array.
    percl (float): The percentile for the lower edge of the stretch (or both edges if ``percu`` is None).
    percu (float): The percentile for the upper edge of the stretch (or None to use ``percl`` for both).
    is_mask (bool): Set to ``True`` if the image is a mask, i.e. all values are either zero or one.
    figsize (tuple): The size of the matplotlib figure in inches.
    cmap (str): Colormap.
    log (bool): If true, use log stretch.
    clip (bool): If true, clip the image.
    show_colorbar (bool): Whether show colorbar or not.
    show_ticks (bool): Whether show ticks or not.
    fig (``matplotlib.pyplot.figure`` object): The figure object.
    ax (``matplotlib.pyplot.axes`` object): The axes object.
    input_ratio (float): The ratio of the input image size to the output image size.

    Returns
    -------
    fig (``matplotlib.pyplot.figure`` object): The figure object.
    """
    if percu is None:
        percu = percl
        percl = 100 - percl

    if (fig is None and ax is not None) or (fig is not None and ax is None):
        raise ValueError(
            'Must provide both "fig" and "ax" ' "if you provide one of them"
        )
    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # To preserve details we should *really* downsample correctly and
    # not rely on matplotlib to do it correctly for us (it won't).

    # So, calculate the size of the figure in pixels, block_reduce to
    # roughly that,and display the block reduced image.

    # Thanks, https://stackoverflow.com/questions/29702424/how-to-get-matplotlib-figure-size
    fig_size_pix = fig.get_size_inches() * fig.dpi

    ratio = (image.shape // fig_size_pix).max()

    if ratio < 1:
        ratio = 1

    ratio = input_ratio or ratio

    # Of course, now that we have downsampled, the axis limits are changed to
    # match the smaller image size. Setting the extent will do the trick to
    # change the axis display back to showing the actual extent of the image.
    extent = [0, image.shape[1], 0, image.shape[0]]

    if log:
        stretch = aviz.LogStretch()
    else:
        stretch = aviz.LinearStretch()

    if vmin is not None and vmax is not None:
        interval = aviz.ManualInterval(vmin, vmax)
    else:
        interval = aviz.AsymmetricPercentileInterval(percl, percu)
    norm = aviz.ImageNormalize(image, interval=interval, stretch=stretch, clip=clip)

    if scale_args is not None:
        scale_args = scale_args
    else:
        if is_mask:
            # The image is a mask in which pixels should be zero or one.
            # block_reduce may have changed some of the values, so reset here.
            image = image > 0
            # Set the image scale limits appropriately.
            scale_args = dict(vmin=0, vmax=1)
        else:
            scale_args = dict(norm=norm)

    im = ax.imshow(
        image,
        origin="lower",
        cmap=cmap,
        extent=extent,
        aspect="equal",
        interpolation=None,
        **scale_args,
    )

    if show_colorbar:
        divider = make_axes_locatable(ax)
        ax_cbar = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, ax=ax, cax=ax_cbar)
        plt.setp(plt.getp(cbar.ax.axes, "xticklabels"), fontsize=11)
        plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), fontsize=11)

    if not show_ticks:
        ax.tick_params(axis="both", which="both", length=0)

    if ax is None:
        return fig
    return fig, ax, scale_args