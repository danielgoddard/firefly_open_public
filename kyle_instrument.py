"""

Provides a set of functions to handle instrumental effects.

:func:`log_rebin` has been pulled from
:mod:`mangadap.contrib.ppxf_util.py` and modified.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/instrument.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import warnings
    import numpy
    from scipy import integrate
    from scipy.interpolate import InterpolatedUnivariateSpline
    from mangadap.util.constants import constants
    from mangadap.util.misc import where_not
    import astropy.constants

*Revision history*:
    | **27 May 2015**: Original implementation by K. Westfall (KBW)
        based on downgrader_MANGA.f provided by D. Thomas, O. Steele, D.
        Wilkinson, D. Goddard.

"""

# from __future__ import division
# from __future__ import print_function
# from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import warnings
import numpy
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from kyle_constants import constants
from kyle_misc import where_not
import astropy.constants
from IPython.core.debugger import Tracer

def spectrum_velocity_scale(wave, log10=False):
    """
    Determine the velocity sampling of an input wavelength coordinate
    vector.  The wavelength vector is assumed to be logarithmically
    sampled, but its units are in angstroms.

    Args: 
        wave (numpy.ndarray): Wavelength coordinates of each spectral
            channel in angstroms.  It is expected that the spectrum has
            been sampled geometrically.
        log10 (bool): (Optional) Input spectrum has been sample
            geometrically using the base 10 logarithm, instead of the
            natural logarithm.

    Returns:
        float : Velocity scale of the spectrum in km/s.

    """

    dl_over_l = (numpy.log10(wave[1])-numpy.log10(wave[0]))*numpy.log(10.0) if log10 else \
                numpy.log(wave[1])-numpy.log(wave[0])
                
    return astropy.constants.c.to('km/s').value*dl_over_l


class convolution_integral_element:
    """
    Support class for variable sigma convolution.  See
    :func:`convolution_variable_sigma`.

    Args:
        y (numpy.ndarray): Vector to convolve
        sigma (numpy.ndarray): Coordinate-dependent standard deviation of the
            Gaussian kernel
        ye (numpy.ndarray): (Optional) Error in the vector to convolve

    Raises:
        Exception: Raised if *y* is not a 1D vector, or if the shape of
            *y* and *sigma* (and *ye* if provided) are different.

    Attributes:
        x (numpy.ndarray): Pixel coordinate vector
        y (numpy.ndarray): Vector to convolve
        ye (numpy.ndarray): Error in the vector to convolve
        sigma (numpy.ndarray): Coordinate-dependent standard deviation of the
            Gaussian kernel
        norm (numpy.ndarray): Gaussian normalization; calculated once for
            efficiency

    """
    def __init__(self, y, sigma, ye=None):
        if len(y.shape) != 1:
            raise Exception('y must be a 1D array!')
        if y.shape != sigma.shape:
            raise Exception('y and sigma must have the same shape!')
        if ye is not None and ye.shape != y.shape:
            raise Exception('y and ye must have the same shape!')
        self.x = numpy.arange(sigma.size, dtype=numpy.float64)
        self.y = y
        self.ye = ye
        self.sigma = sigma
        self.norm = numpy.sqrt(2.0*numpy.pi)*self.sigma

    
    def _get_kernel(self, xc):
        """Calculate the kernel vector when centered at *xc*."""
        return numpy.exp(-0.5*numpy.square((self.x-xc)/self.sigma))/self.norm

    def __call__(self, xc):
        """
        Calculates the weighted mean of :attr:`y`, where the weights are
        defined by a Gaussian with standard deviation :attr:`sigma` and
        centered at xc.

        Args:
            xc (float): Center for the Gaussian weighting function

        Returns:
            float: The weighted mean of :attr:`y`

        """
        kernel = self._get_kernel(xc)
        return numpy.sum(self.y*kernel) / numpy.sum(kernel)


    def error(self, xc):
        """
        Calculates the error in the weighted mean of :attr:`y` using
        nominal error propagation.  The weights are defined by a
        Gaussian with standard deviation :attr:`sigma` and centered at
        xc.

        Args:
            xc (float): Center for the Gaussian weighting function

        Returns:
            float: The error in the weighted mean of :attr:`y`

        """
        kernel = self._get_kernel(xc)
        return numpy.sqrt(numpy.sum(numpy.square(self.ye*kernel)) / numpy.sum(kernel))


def convolution_variable_sigma(y, sigma, ye=None):
    r"""
    Convolve a discretely sampled function :math:`y(x)` with a Gaussian
    kernel, :math:`g`, where the standard deviation of the kernel is a
    function of :math:`x`, :math:`\sigma(x)`.  Nominal calculations can
    be performed to propagate the error in the result; these
    calculations **do not** include the covariance between the pixels,
    which will mean that the calculations likely have significant error!

    The convolution is defined as:

    .. math::

        (y\ast g)(x) &= \int_{-\infty}^{\infty} y(X)\ g(\sigma,x-X)\ dX \\
                     &= \int_{-\infty}^{\infty} \frac{y(X)}{\sqrt{2\pi}\
                        \sigma(X)}\ \exp\left(-\frac{(x-X)^2}{2\
                        \sigma(X)^2}\right) dX .

    To minimize edge effects and account for the censoring of the data
    (finite range in :math:`x`), the convolution is actually calculated
    as a definite integral and normalized as follows:

    .. math::

        (y\ast g)(x) \sim\frac{
        \int_{x-n_\sigma*\sigma(x)}^{x+n_\sigma*\sigma(x)} y(X)\
        g(\sigma,x-X)\ dX}{
        \int_{x-n_\sigma*\sigma(x)}^{x+n_\sigma*\sigma(x)}
        g(\sigma,x-X)\ dX} .

    The above is identical to getting the weighted mean of :math:`y` at
    each position :math:`x`, where the weights are defined by a Gaussian
    kernel centered at :math:`x` with a variable dispersion.

    Use of this function requires:
        - *y* and *sigma* must be 1D vectors
        - *y* and *sigma* must be uniformly sampled on the same grid
        - *sigma* must be in pixel units.

    Args:
        y (numpy.ndarray): A uniformly sampled function to convolve.
        sigma (numpy.ndarray): The standard deviation of the Gaussian
            kernel sampled at the same positions as *y*.  The units of
            *sigma* **must** be in pixels.
        ye (numpy.ndarray): (Optional) Errors in the function
            :math:`y(x)`.

    Returns:
        numpy.ndarray: Arrays with the convolved function :math:`(y\ast
            g)(x)` sampled at the same positions as the input :math:`x`
            vector and its error.  The second array will be returned as
            None if the error vector is not provided.

    Raises:
        Exception: Raised if trying to calculate the errors because they
        haven't been implemented yet.

    """
    kernel = convolution_integral_element(y,sigma,ye=ye)

    #numpy.exp(-0.5*numpy.square((self.x-xc)/self.sigma))/self.norm
    conv = numpy.array(list(map(kernel, kernel.x)))
    print numpy.shape(y)
    # print sys.version
    # print "hi"
    # Tracer()()

    if ye is None:
        return conv

    out = numpy.array(list(map(kernel.error, kernel.x)))
    return conv, out


class spectral_resolution:
    r"""
    
    Container class for the resolution, :math:`R =
    \lambda/\Delta\lambda`, of a spectrum.  The primary functionality is
    to determine the parameters necessary to match the resolution of one
    spectrum to another.  It can also be used as a function to
    interpolate the spectral resolution at a given wavelenth.

    Args:
        wave (numpy.ndarray): A 1D vector with the wavelength in
            angstroms.  The sampling may be either in linear steps of
            wavelength or :math:`\log_{10}` steps.
        sres (numpy.ndarray): A 1D vector with the spectral resolution,
            :math:`R`, sampled at the positions of the provided
            wavelength vector.
        log10 (bool): (Optional) Flag that the spectrum has been binned
            logarithmically (base 10) in wavelength
        interp_ext (int or str): (Optional) The value to pass as *ext*
            to the interpolator, which defines its behavior when
            attempting to sample the spectral resolution beyond where it
            is defined.  See
            `scipy.interpolate.InterpolatedUnivariateSpline`_.  Default
            is to extrapolate.

    Raises:
        Exception: Raised if *wave* is not a 1D vector or if *wave* and
            *sres* do not have the same shape.

    Attributes:
        interpolator
            (`scipy.interpolate.InterpolatedUnivariateSpline`_): An
            object used to interpolate the spectral resolution at any
            given wavelength.  The interpolation is hard-wired to be
            **linear** and its extrapolation behavior is defined by
            *interp_ext*.  The wavelength and resolution vectors are
            held by this object for later reference if needed.
        log10 (bool): Flag that the spectrum has been binned
            logarithmically (base 10) in wavelength
        cnst (:class:`mangadap.util.constants`): Used to define the
            conversion factor between a Gaussian sigma and FWHM
        c (float): The speed of light; provided by `astropy.constants`_.
        dv (float): The velocity step per pixel in km/s.  Defined using
            :func:`spectrum_velocity_scale` if :attr:`log10` is True;
            otherwise set to None.
        dw (float): The wavelength step per pixel in angstroms.  Defined
            as the wavelength step between the first two pixels if
            :attr:`log10` is False; otherwise set to None.
        min_sig (float): Minimum standard deviation allowed for the
            kernel used to match two spectral resolutions.  See
            :func:`GaussianKernelDifference`.
        sig_pd (numpy.ndarray): The standard deviation in pixels
            required to match the spectral resolution of this object to
            the resolution defined by a different spectral_resolution
            object.  See :func:`GaussianKernelDifference`.
        sig_mask (numpy.ndarray): A *uint* vector used to identify
            measurements of :attr:`sig_pd` that should **not** be used
            to match the spectral resolution of this object to the
            resolution defined by a different spectral_resolution
            object.  See :func:`GaussianKernelDifference`.
        sig_vo (float): A constant offset of the kernal standard
            deviation **in km/s** that has been applied to
            :attr:`sig_pd`.  See :func:`GaussianKernelDifference`.
        
    .. todo::

        - Allow it to determine if the binning is linear or geometric,
          then use the *log10* keyword to distinguish between natural
          log and :math:`log_{10}` binning.
        - Allow for more than one type of line-spread function?

    .. warning::

        The default behavior of the interpolator is to extrapolate the
        input spectral resolution vector when trying to sample from
        regions beyond where it is sampled.  Use *interp_ext* change
        this; see `scipy.interpolate.InterpolatedUnivariateSpline`_.

    .. _scipy.interpolate.InterpolatedUnivariateSpline: http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.InterpolatedUnivariateSpline.html
    .. _astropy.constants: http://docs.astropy.org/en/stable/constants/index.html

    """
    def __init__(self, wave, sres, log10=False, interp_ext='extrapolate'):
        # Check the sizes
        if len(wave.shape) != 1:
            raise Exception('wave must be a 1D array!')
        if wave.shape != sres.shape:
            raise Exception('wave and sres must have the same shape!')

        # k=1; always use linear interpolation
        self.interpolator = InterpolatedUnivariateSpline(wave, sres, k=1, ext=interp_ext)
        self.log10 = log10
        self.cnst = constants()
        self.c = astropy.constants.c.to('km/s').value

        self.dv = spectrum_velocity_scale(wave, log10=True) if log10 else None
        self.dw = None if log10 else wave[1] - wave[0]

        # No resolution matching calculated yet
        self.min_sig = None
        self.sig_pd = None
        self.sig_mask = None
        self.sig_vo = None


    def __call__(self, w):
        """Interpolate the spectral resolution at wavelength *w*."""
        return self.interpolator(w)


    def _finalize_GaussianKernelDifference(self, sig2_pd):
        r"""
        Given the calculated :math:`\sigma^2_{p,d}`, calculate and save
        the attributes :attr:`sig_pd` and :attr:`sig_mask`.  See
        :func:`GaussianKernelDifference`.
        """
        indx = numpy.where(numpy.isclose(sig2_pd, 0.0))
        nindx = where_not(indx, sig2_pd.size)
        self.sig_pd = numpy.copy(sig2_pd)
        self.sig_pd[nindx] = sig2_pd[nindx]/numpy.sqrt(numpy.absolute(sig2_pd[nindx]))
        self.sig_pd[indx] = 0.0
        self.sig_mask = numpy.array(self.sig_pd < -self.min_sig).astype(numpy.uint)


    def _convert_vd2pd(self, sig2_vd):
        r"""
        Convert from :math:`\sigma^2_{v,d}` to :math:`\sigma^2_{p,d}`.
        See :func:`GaussianKernelDifference`.
        """
        return sig2_vd / numpy.square(self.dv) if self.log10 else \
               sig2_vd / numpy.square(self.c*self.dw/self.wave())


    def _convert_pd2vd(self, sig2_pd):
        r"""
        Convert from :math:`\sigma^2_{p,d}` to :math:`\sigma^2_{v,d}`.
        See :func:`GaussianKernelDifference`.
        """
        return sig2_pd * numpy.square(self.dv) if self.log10 else \
               sig2_pd * numpy.square(self.c*self.dw/self.wave())


    def wave(self):
        """
        Return the wavelength vector; held by :attr:`interpolator`.
        """
        return self.interpolator._data[0]


    def sres(self):
        """
        Return the resolution vector; held by :attr:`interpolator`.
        """
        return self.interpolator._data[1]


    def match(self, new_sres, no_offset=True, min_sig_pix=0.0):
        """
        Currently only an alias for :func:`GaussianKernelDifference`.
        """
        self.GaussianKernelDifference(new_sres, no_offset=no_offset, min_sig_pix=min_sig_pix)


    def GaussianKernelDifference(self, new_sres, no_offset=True, min_sig_pix=0.0):
        r"""
        Determine the parameters of a Gaussian kernel required to
        convert the resolution of this object to the resolution of a
        different the :class:`spectral_resolution` object, *new_sres*.

        The spectral resolution is defined as :math:`R =
        \lambda/\Delta\lambda`, where :math:`\Delta\lambda` is the FWHM
        of the spectral resolution element.  The standard deviation of
        the resolution element in angstroms is then
   
        .. math::
    
            \sigma_\lambda = \frac{\lambda}{f R}, \ \ {\rm where} \ \  f
            = \frac{{\rm FWHM_\lambda}}{\sigma_\lambda}.

        Assuming a Gaussian (in angstroms) line-spread function:

        .. math::

            \sigma^2_{\lambda,2} = \sigma^2_{\lambda,1} +
            \sigma^2_{\lambda,d}

        such that

        .. math::

            \sigma^2_{\lambda,d} = \left(\frac{\lambda}{f}\right)^2
            (R^{-2}_2 - R^{-2}_1)

        is the defining parameter of the Gaussian kernel needed to take
        a spectrum of resolution :math:`R_1` to one with a resolution of
        :math:`R_2`.

        For input to :func:`convolution_variable_sigma`, the
        *wavelength-dependent* parameter of the Gaussian kernel is
        converted to pixels.  This function allows for the input spectra
        to be linearly sampled in angstroms or log10(angstroms).  For
        the former (*log10=False*), 

        .. math::

            \sigma^2_{p,d} = \left(\frac{\lambda}{f\
            \delta\lambda}\right)^2 (R^{-2}_2 - R^{-2}_1)

        where :math:`\delta\lambda` is the size of the pixel in
        angstroms.  If the units are log10(angstrom) (*log10=True*), we
        approximate the velocity width of each pixel to be :math:`\delta
        v = c \ln(10.0) (\log\lambda[1]-\log\lambda[0])`, such that
    
        .. math::

            \sigma^2_{p,d} &= \left(\frac{c}{ \delta v \lambda}\right)^2
            \sigma^2_{\lambda,d} \\ &= \left(\frac{c}{ f\ \delta
            v}\right)^2 (R^{-2}_2 - R^{-2}_1)\ ;

        :math:`c` is the speed of light in km/s.

        The nominal use of this algorithm assumes :math:`R_1 \geq R_2`.
        However, in practice, :func:`convolution_variable_sigma` only
        uses a Gaussian kernel up to some minimum value of
        :math:`\epsilon_\sigma`; below this, the kernel is assumed to be
        a Delta function.  Therefore, as long as
    
        .. math::
    
            \sigma_{p,d} \equiv \sigma^2_{p,d}/\sqrt{|\sigma^2_{p,d}|}
            \geq -\epsilon_\sigma\ ,
        
        the behavior of :func:`convolution_variable_sigma` should not be
        affected.
    
        Even so, there may be spectral regions that do not have
        :math:`\sigma_{p,d} \geq -\epsilon_\sigma`; for such spectral
        regions there are three choices:

            (**Option 1**) trim the spectral range to only those
            spectral regions where the existing resolution is better
            than the target resolution,
        
            (**Option 2**) match the existing resolution to the target
            resolution up to some constant offset that must be accounted
            for in subsequent analyses, or

            (**Option 3**) allow for a wavelength dependent difference
            in the spectral resolution that must be accounted for in
            subsequent analyses.

        The choice of either Option 1 or 2 is selected by setting
        *no_offset* to, respectively, True or False; Option 1 is the
        default behavior.  Currently, Option 3 is not allowed.

        For Option 1, pixels with :math:`\sigma_{p,d} <
        -\epsilon_\sigma` are masked (*sigma_mask = 1*); however, the
        returned values of :math:`\sigma_{p,d}` are left unchanged.

        For Option 2, we define

        .. math::

            \sigma^2_{v,o} = -{\rm min}(\sigma^2_{v,d}) - {\rm
            max}(\epsilon_\sigma \delta v)^2

        where :math:`\delta v` is constant for the logarithmically
        binned spectrum and is wavelength dependent for the linearly
        binned spectra; in the latter case, the velocity step is
        determined for each pixel::

            _wave = self.wave()
            dv = self.c * (2.0*(_wave[1:] - _wave[0:-1]) / (_wave[1:] + _wave[0:-1]))

        If :math:`\sigma^2_{v,o} > 0.0`, it must be that :math:`{\rm
        min}(\sigma^2_{v,d}) < -{\rm max}(\epsilon_\sigma \delta v)^2`,
        such that an offset should be applied.  In that case, the
        returned kernel parameters are

        .. math::

            \sigma^\prime_{v,d} = \sqrt{\sigma^2_{v,d} +
            \sigma^2_{v,o}}\ .

        with the units converted to pixels using the equations above, no
        pixels are masked, and :math:`\sqrt{\sigma^2_{v,o}}` is returned
        for the offset.  Otherwise, the offset is set to 0.

        .. todo::

            Allow to check cases when the convolution kernel is
            indpendent of wavelength such that the convolution can be
            sped up by performing the convolution using an FFT.  For
            example, in the case where the spectrum is logarithmically
            binned and both :math:`R_1` and :math:`R_2` are
            *independent* of wavelength, the convolution kernel is
            independent of wavelength.

        Args:
            new_sres (:class:`spectral_resolution`): Spectral resolution
                to match to.

            no_offset (bool): (Optional) Force :math:`\sigma^2_{v,o} =
                0` by masking regions with :math:`\sigma_{p,d} <
                -\epsilon_\sigma`; i.e., the value of this arguments
                selects Option 1 (True) or Option 2 (False).

            min_sig_pix (float): (Optional) Minimum value of the
                standard deviation allowed before assuming the kernel is
                a Delta function.

        """
        # Save the minimum pixel sigma to allow
        self.min_sig = min_sig_pix

        # Interpolate the new spectral resolution vector at the wavelengths
        # of the input spectral resolution
        _wave = self.wave()
        _sres = self.sres()
        interp_sres = new_sres(_wave)

        # Determine the variance (in angstroms) of Gaussian needed to match
        # input resolution to the new values
        sig2_wd = numpy.square(_wave/self.cnst.sig2fwhm) \
                  * (1.0/numpy.square(interp_sres) - 1.0/numpy.square(_sres))
        # Convert to km/s
        sig2_vd = numpy.square(self.c/_wave) * sig2_wd

        # Option 1:
        if no_offset:
            # Convert the variance to pixel coordinates
            sig2_pd = sig2_vd / numpy.square(self.dv) if self.log10 else \
                      sig2_wd / numpy.square(self.dw)
            # No offset applied
            self.sig_vo = 0.0

        # Option 2:
        else:
            # Calculate the velocity step of each pixel
            dv = self.c * (2.0*(_wave[1:] - _wave[0:-1]) / (_wave[1:] + _wave[0:-1]))
            # Get the needed *velocity* offset (this is the square)
            self.sig_vo = - numpy.amin(sig2_vd) - numpy.square(self.min_sig * numpy.amax(dv))
            # Apply it if it's larger than 0
            if self.sig_vo > 0:
                sig2_vd += self.sig_vo
                self.sig_vo = numpy.sqrt(self.sig_vo)
            else:
                self.sig_vo = 0.0

            # Convert the variance to pixel coordinates
            sig2_pd = self._convert_vd2pd(sig2_vd)

        self._finalize_GaussianKernelDifference(sig2_pd)


    def offset_GaussianKernelDifference(self, offset):
        r"""
        If the properties required to match the resolution of one
        spectrum to another has already been calculated (see
        :func:`GaussianKernelDifference`), this allows for one to apply
        an additional offset.  The additional offset **must** be in km/s
        (not pixels).

        The offset is applied in quadrature; however, the offset can be
        negative such that one can reduce :attr:`sig_pd`.  Once
        converted to km/s, the offset is applied by calculating:

        .. math::
        
            \sigma^{\prime\ 2}_{v,d} = \sigma^{2}_{v,d} +
            \sigma_{off}|\sigma_{off}|\ .

        :attr:`sig_vo` is adjusted in the same way, and the change in
        :math:`\sigma^{\prime\ 2}_{v,d}` is then propagated to
        :attr:`sig_pd` and :attr:`sig_mask`.
        
        Args:
            offset (float): Value of the standard deviation in km/s to
                add in quadrature to a previously calculated
                :attr:`sig_pd`.
        
        Raises:
            Exception: Raised if the kernel properties have not yet been
                defined.

        """
        if None in [self.min_sig, self.sig_pd, self.sig_mask, self.sig_vo]:
            raise Exception('No kernel defined yet.  Run GaussianKernelDifference first.')
        if numpy.isclose(offset,0.0):
            return
        off2 = offset*numpy.absolute(offset)
        sig2_vo = self.sig_vo*numpy.absolute(self.sig_vo) + off2
        self.sig_vo = 0.0 if numpy.isclose(sig2_vo, 0.0) \
                          else sig2_vo/numpy.sqrt(numpy.absolute(sig2_vo))
        sig2_vd = self._convert_pd2vd(self.sig_pd*numpy.absolute(sig_pd)) + off2
        self._finalize_GaussianKernelDifference(self._convert_vd2pd(sig2_vd))


    def adjusted_resolution(self, indx=None):
        r"""

        Return the resolution that should result from applying
        :func:`convolution_variable_sigma` to the spectrum associated
        with this spectral resolution object using :attr:`sigma_pd`.
        I.e., calculate:

        .. math::

            R_2 = \left[ \left(\frac{f}{c}\right)^2 \sigma^2_{v,d} +
            R^{-2}_1\right]^{-1/2}\ . 

        Args:
            indx (tuple): (Optional) Selection tuple used to return a
                subset of the full resolution vector.

        Returns:
            numpy.ndarray: The (full or selected) vector with the
                adjusted resolution.

        .. todo::
            Allow to reset the resolution of this object to the adjusted
            resolution and reset the kernel variables to None.

        """
        if indx is None:
            return 1.0/numpy.sqrt( numpy.square(self.cnst.sig2fwhm/self.c) \
                                   * self._convert_pd2vd(self.sig_pd*numpy.absolute(self.sig_pd)) \
                                   + 1.0/numpy.square(self.sres()) )

        return 1.0/numpy.sqrt( numpy.square(self.cnst.sig2fwhm/self.c) \
                            * (self._convert_pd2vd(self.sig_pd*numpy.absolute(self.sig_pd)))[indx] \
                               + 1.0/numpy.square(self.sres()[indx]) )


def match_spectral_resolution(wave, flux, sres, new_sres_wave, new_sres, ivar=None, mask=None, 
                              min_sig_pix=0.0, no_offset=True, variable_offset=False, log10=False,
                              new_log10=False):
    r"""
    Adjust the existing spectral resolution of to a **lower** resolution
    as best as possible.  The primary functionality is in
    :class:`spectral_resolution`, which determines the Gaussian kernel
    parameters needed to match the resolution, and
    :func:`convolve_variable_sigma`, which actually performs the
    convolution to match the resolution.

    In particular, see
    :func:`spectral_resolution.GaussianKernelDifference` for a
    description of how the kernel parameters are determined and how
    regions where the target resolution is **higher** than the existing
    resolution.  In this case, one of the options is to adopt an offset
    of the resolution (in km/s) that could be corrected for in
    subsequent analysis.  In this case, setting *variable_offset* to
    True allows the offset to be different for all input spectra.  If
    one expects to combine the spectra, the default behavior should be
    used, forcing all the spectra to have a constant offset.

    Args:
        wave (numpy.ndarray): A 1D or 2D (:math:`N_{\rm spec}\times
            N_{\rm pix}`) array with the wavelength in angstroms for a
            set of spectra.  The sampling may be either in linear steps
            of wavelength or :math:`\log_{10}` steps, as set using
            *log10*.
        flux (numpy.ndarray): A 1D or 2D (:math:`N_{\rm spec}\times
            N_{\rm pix}`) array with the flux sampled at the provided
            wavelengths.
        sres (numpy.ndarray): A 1D or 2D (:math:`N_{\rm spec}\times
            N_{\rm pix}`) array with the spectral resolution, :math:`R`,
            at the provided wavelengths.
        new_sres_wave (numpy.ndarray): A 1D vector with the wavelength
            in angstroms at which the new resolution of the input
            spectra has been sampled.  The sampling may be either in
            linear steps of wavelength or :math:`\log_{10}` steps, as
            set using *new_log10*.
        new_sres (numpy.ndarray): A 1D vector with the new resolution
            for the input spectra.
        ivar (numpy.ndarray): (Optional) A 1D or 2D (:math:`N_{\rm
            spec}\times N_{\rm pix}`) array with the inverse variance of
            the flux sampled at the provided wavelengths.  **Currently
            never used and, if provided, the output inverse variances
            are simply a carbon copy of this array.**
        mask (numpy.ndarray): (Optional) A 1D or 2D (:math:`N_{\rm
            spec}\times N_{\rm pix}`) array with a *uint* mask for the
            flux sampled at the provided wavelengths.
        no_offset (bool): (Optional) Force :math:`\sigma^2_{v,o} = 0` by
            masking regions with :math:`\sigma_{p,d} <
            -\epsilon_\sigma`; i.e., the value of this arguments selects
            Option 1 (True) or Option 2 (False).  See
            :func:`spectral_resolution.GaussianKernelDifference`.
        min_sig_pix (float): (Optional) Minimum value of the standard
            deviation in pixels allowed before assuming the kernel is a
            Delta function.
        variable_offset (bool): (Optional) Flag to allow the offset
            applied to each spectrum (when the input contains more than
            one spectraum) to be tailored to each spectrum.  Otherwise
            (*variable_offset=False*) the offset is forced to be the
            same for all spectra.
        log10 (bool): (Optional) Flag that the spectrum has been binned
            logarithmically (base 10) in wavelength
        new_log10 (bool): (Optional) Flag that the coordinates of the
            new spectral resolution are  spectrum as been binned
            logarithmically (base 10) in wavelength.

    Returns: 

        numpy.ndarray : Four or Five arrays are returned:

            - A 1D or 2D (:math:`N_{\rm spec}\times N_{\rm pix}`) array
              with the resolution-matched flux sampled at the input
              wavelengths.
            - A 1D or 2D (:math:`N_{\rm spec}\times N_{\rm pix}`) array
              with the spectral resolution, :math:`R`, of the
              resolution-matched spectra at the provided wavelengths.
            - A 1D vector with any constant offset in resolution **in
              km/s** between the targetted value and the result.  See
              :func:`spectral_resolution.GaussianKernelDifference`.
            - A 1D or 2D (:math:`N_{\rm spec}\times N_{\rm pix}`) array
              with a *uint* mask for the resolution-matched flux sampled
              at the input wavelengths.  This is returned regardless of
              whether an input mask was provided.  Any pixel that had a
              resolution that was lower than the target resolution (up
              to some tolerance defined by *min_sig_pix*) is returned as
              masked.
            - (Optional) A 1D or 2D (:math:`N_{\rm spec}\times N_{\rm
              pix}`) array with the inverse variance of the
              resolution-matched flux sampled at the input wavelengths.
              **Only output if *ivar* is provided, and in that case this
              is just a carbon copy of the input array.**

    Raises:
        Exception: Raised if:

            - the input *wave* array is 2D and the *sres* array is not;
              a 1D wavelength array is allowed for a 2D *sres* array but
              not vice versa

            - the number of spectral pixels in *wave*, *flux*, and
              *sres* is not the same

            - the shape of the *flux*, *mask* (if provided), and *ivar*
              (if provided) are not the same

            - the shape of the *new_sres_wave* and *new_sres* arrays
              are not the same and/or not 1D

    .. todo::

        - Add interp_ext != 'extrapolate' option?
        - Better way to use warnings?

    """
    # Check the dimensionality of wave and sres
    wave_matrix = len(wave.shape) == 2
    sres_matrix = len(sres.shape) == 2
    if wave_matrix and not sres_matrix:
        raise Exception('If input wavelength array is 2D, the spectral resolution array must' \
                        ' also be 2D')

    # Check the shapes
    if (wave_matrix == sres_matrix and wave.shape != sres.shape) or \
       (not wave_matrix and sres_matrix and wave.shape[0] != sres.shape[1]):
        raise Exception('Input spectral resolution and coordinate arrays must have the same' \
                        ' number of spectral channels!')
    if (wave_matrix and wave.shape != flux.shape) or \
       (not wave_matrix and len(flux.shape) == 2 and wave.shape[0] != flux.shape[1]) or \
       (not wave_matrix and len(flux.shape) == 1 and wave.shape != flux.shape):
        raise Exception('Input flux and coordinate arrays must have the same number of' \
                        ' spectral channels!')
    if (mask is not None and mask.shape != flux.shape):
        raise Exception('Input flux and mask arrays must have the same shape!')
    if (ivar is not None and ivar.shape != flux.shape):
        raise Exception('Input flux and ivar arrays must have the same shape!')
        
    if len(new_sres_wave.shape) != 1 or len(new_sres.shape) != 1:
        raise Exception('New spectral resolution and coordinate arrays must be 1D!')
    if new_sres_wave.shape != new_sres.shape:
        raise Exception('New spectral resolution and coordinate arrays must have the same shape!')

    # Raise a warning if the new_sres vector will have to be
    # extrapolated for the input wavelengths
    if numpy.amin(wave) < new_sres_wave[0] or numpy.amax(wave) > new_sres_wave[-1]:
        warnings.warn('Mapping to the new spectral resolution will require extrapolating the' \
                      ' provided input vectors!')

    # Initialize some variables
    nspec = 1 if len(sres.shape) == 1 else sres.shape[0]
    sigma_offset = numpy.zeros(nspec, dtype=numpy.float64)
    new_res = spectral_resolution(new_sres_wave, new_sres, log10=new_log10)
    res = []

    # Get the kernel parameters necessary to match all spectra to the
    # new resolution
    if nspec == 1:
        res = [spectral_resolution(wave, sres, log10=log10)]
        res[0].match(new_res, no_offset=no_offset, min_sig_pix=min_sig_pix)
        sigma_offset[0] = res[0].sig_vo
    else:
        for i in range(0,nspec):
            _wave = wave[i,:].ravel() if wave_matrix else wave
            _sres = sres[i,:].ravel() if sres_matrix else sres
            res = res + [spectral_resolution(_wave, _sres, log10=log10)]
            res[i].match(new_res, no_offset=no_offset, min_sig_pix=min_sig_pix)
            sigma_offset[i] = res[i].sig_vo

    # Force all the offsets to be the same, if requested
    if not no_offset and not variable_offset:
        common_offset = numpy.max(sigma_offset)
        offset_diff = numpy.sqrt( numpy.square(common_offset) - numpy.square(sigma_offset))
        for r, o in zip(res,offset_diff):
            r.offset_GaussianKernelDifference(o)

    # Perform the convolutions
    out_flux = numpy.copy(flux)
    out_sres = numpy.copy(sres)
    if mask is None:
        mask = numpy.zeros(flux.shape, dtype=numpy.uint)
    out_mask = numpy.copy(mask)

    if nspec == 1:
        indx = numpy.where(res[0].sig_pd > min_sig_pix)
        out_flux[indx] = convolution_variable_sigma(flux[indx], res[0].sig_pd[indx])
        out_sres[indx] = res[0].adjusted_resolution(indx=indx)
        out_mask = numpy.array((res[0].sig_mask == 1) | (mask == 1)).astype(numpy.uint)
    else:
        for i in range(0,nspec):
            print('Matching resolution: {0}/{1}'.format(i+1,nspec))
            indx = numpy.where(res[i].sig_pd > min_sig_pix)
            out_flux[i,indx] = convolution_variable_sigma(flux[i,indx].ravel(), res[i].sig_pd[indx])
            out_sres[i,indx] = res[i].adjusted_resolution(indx=indx)
            out_mask[i,:] = numpy.array((res[i].sig_mask == 1) \
                                        | (mask[i,:] == 1)).astype(numpy.uint)

    # TODO: Add this functionality from the IDL version?
    #
    # Finally, the code masks a number of pixels at the beginning and
    # end of the spectra to remove regions affected by errors in the
    # convolution due to the censoring of the data.  The number of
    # pixels is the FWHM of the largest Gaussian applied in the
    # convolution: ceil(sig2fwhm*max(diff_sig_w)/dw).  This is currently
    # hard-wired and should be tested.

    if ivar is not None:
        out_ivar = numpy.copy(ivar)
        # Handle the inverse variances then return the answers
        return out_flux, out_sres, sigma_offset, out_mask, out_ivar

    return out_flux, out_sres, sigma_offset, out_mask
    



def log_rebin(lamRange, spec, oversample=None, velscale=None, flux=False, log10=False,
              newRange=None, wave_in_ang=False, unobs=0.0):
    """
    .. note::
     
        Copyright (C) 2001-2014, Michele Cappellari
        E-mail: cappellari_at_astro.ox.ac.uk
     
        This software is provided as is without any warranty whatsoever.
        Permission to use, for non-commercial purposes is granted.
        Permission to modify for personal or internal use is granted,
        provided this copyright and disclaimer are included unchanged at
        the beginning of the file. All other rights are reserved.

    Logarithmically rebin a spectrum, while rigorously conserving the
    flux.  Basically the photons in the spectrum are simply
    ridistributed according to a new grid of pixels, with non-uniform
    size in the spectral direction.
    
    This routine makes the `standard' zero-order assumption that the
    spectrum is *constant* within each pixels. It is possible to perform
    log-rebinning by assuming the spectrum is represented by a
    piece-wise polynomial of higer degree, while still obtaining a
    uniquely defined linear problem, but this reduces to a deconvolution
    and amplifies noise.

    .. warning::

        This assumption can be poor for sharp features in the spectrum.
        Beware if resampling spectra with strong, marginally sampled
        features!
    
    This same routine can be used to compute approximate errors of the
    log-rebinned spectrum. To do this type the command
    
    >>> err2New, logLam, velscale = log_rebin(lamRange, numpy.square(err))
    
    and the desired errors will be given by numpy.sqrt(err2New).
    
    .. warning::
    
        This rebinning of the error-spectrum is very *approximate* as it
        does not consider the correlation introduced by the rebinning!

    Args:

        lamRange (numpy.ndarray): two elements vector containing the
            central wavelength of the first and last pixels in the
            spectrum, which is assumed to have constant wavelength
            scale! E.g. from the values in the standard FITS keywords:
            LAMRANGE = CRVAL1 + [0,CDELT1*(NAXIS1-1)].  It must be
            LAMRANGE[0] < LAMRANGE[1].
        spec (numpy.ndarray): Input spectrum.
        oversample (int): (Optional) Oversampling can be done, not to
            loose spectral resolution, especally for extended wavelength
            ranges and to avoid aliasing.  Default is to provide the
            same number of output pixels as input.
        velscale (float): Velocity scale in km/s per pixels. If this
            variable is not defined, then it will contain in output the
            velocity scale.  If this variable is defined by the user it
            will be used to set the output number of pixels and
            wavelength scale.
        flux (bool): (Optional) Set this keyword to preserve total flux.
            In this case the log rebinning changes the pixels flux in
            proportion to their dLam so the following command will show
            large differences beween the spectral shape before and after
            :func:`log_rebin`::
     
                # Plot log-rebinned spectrum
                pyplot.plot(exp(logLam), specNew)
                pyplot.plot(numpy.arange(lamRange[0],lamRange[1],spec.size), spec, 'g')
                pyplot.show()
     
            By default, when this keyword is *not* set, the above two
            lines produce two spectra that almost perfectly overlap each
            other.
        log10 (bool): (Optional) Flag that the spectrum should be binned
            in units of base-10 log wavelength, instead of natural log
        newRange (numpy.ndarray): (Optional) Force the spectrum to be
            sampled to a new spectral range (lamRange is the *existing*
            spectral range).
        wave_in_ang (bool): (Optional) Return the wavelength coordinates
            in angstroms, not log(angstroms)
        unobs (float): (Optional) Default value for unobserved spectral
            regions.

    Returns:
        numpy.ndarray, float : Returns three variables: logarithmically
            rebinned spectrum, the log of the wavelength at the
            geometric center of each pixel, and the velocity scale of
            each pixel in km/s.
        
    Raises:
        ValueError: Raised if the input spectrum is not a
            one-dimensional numpy.ndarray.
        
    *Modification History*:
        | **V1.0.0**: Using interpolation. Michele Cappellari, Leiden,
            22 October 2001
        | **V2.0.0**: Analytic flux conservation. MC, Potsdam, 15 June
            2003
        | **V2.1.0**: Allow a velocity scale to be specified by the
            user.  MC, Leiden, 2 August 2003
        | **V2.2.0**: Output the optional logarithmically spaced
            wavelength at the geometric mean of the wavelength at the
            border of each pixel.  Thanks to Jesus Falcon-Barroso. MC,
            Leiden, 5 November 2003
        | **V2.2.1**: Verify that lamRange[0] < lamRange[1].  MC,
            Vicenza, 29 December 2004
        | **V2.2.2**: Modified the documentation after feedback from
            James Price.  MC, Oxford, 21 October 2010
        | **V2.3.0**: By default now preserve the shape of the spectrum,
            not the total flux. This seems what most users expect from
            the procedure.  Set the keyword /FLUX to preserve flux like
            in previous version.  MC, Oxford, 30 November 2011
        | **V3.0.0**: Translated from IDL into Python. MC, Santiago, 23
            November 2013
        | **V3.1.0**: Fully vectorized log_rebin. Typical speed up by
            two orders of magnitude.  MC, Oxford, 4 March 2014
        | **05 Jun 2015**: (K. Westfall, KBW) Pulled from ppxf_util.py.
            Conform to mangadap documentation standard.  Transcribe
            edits made to IDL version that provides for the log10 and
            newRange arguments.  Add option to return wavelength in
            angstroms, not log(angstroms).  Break out determination of
            input and output spectrum pixel coordinates to a new
            function, :func:`log_rebin_pix`.  Added default value for
            unobserved pixels.  Default behavior unchanged.

    .. todo::

        - Allow to resample an already geometrically binned spectrum
    
    """
    lamRange = numpy.asarray(lamRange)

    if type(spec) != numpy.ndarray:
        raise ValueError('Input spectrum must be a numpy.ndarray')
    s = spec.shape
    if len(s) != 1:
        raise ValueError('input spectrum must be a vector')
    n = s[0]

    # This is broken out into a separate procedure so that it can be
    # called to determine the size of the rebinned spectra without
    # actually doing the rebinning
    dLam, m, logscale, velscale = \
        log_rebin_pix(lamRange, n, oversample=oversample, velscale=velscale, log10=log10,
                      newRange=newRange)

    # Get the sampling of the existing spectrum
    lim = lamRange/dLam + [-0.5, 0.5]           # All in units of dLam
    borders = numpy.linspace(*lim, num=n+1)     # Linearly sampled pixels

    # Set limits to a new wavelength range
    if newRange is not None:
        lim = numpy.asarray(newRange)/dLam + [-0.5, 0.5]

    # Set the limits to the (base-10 or natural) log of the wavelength
    logLim = numpy.log(lim) if not log10 else numpy.log10(lim)
    logLim[1] = logLim[0] + m*logscale      # Set last wavelength, based on integer # of pixels

    # Geometrically spaced pixel borders for the new spectrum
    newBorders = numpy.power(10., numpy.linspace(*logLim, num=m+1)) if log10 else \
                 numpy.exp(numpy.linspace(*logLim, num=m+1))

    # Get the new spectrum by performing an analytic integral
    k = (newBorders - borders[0]).clip(0, n-1).astype(int)
    specNew = numpy.add.reduceat(spec, k)[:-1]
    specNew *= numpy.diff(k) > 0                # fix for design flaw of reduceat()
    specNew += numpy.diff((newBorders - borders[k])*spec[k])

    # Don't conserve the flux
    if not flux:
        specNew /= numpy.diff(newBorders)

    # Output log(wavelength): log of geometric mean
    LamNew = numpy.sqrt(newBorders[1:]*newBorders[:-1])*dLam

    # Set values for unobserved regions
    if newRange is not None and (newRange[0] < lamRange[0] or newRange[1] > lamRange[1]):
            specNew[ (LamNew < lamRange[0]) | (LamNew > lamRange[1]) ] = unobs

    # Return log(wavelength), if requested
    if not wave_in_ang:
        LamNew = numpy.log10(LamNew) if log10 else numpy.log(LamNew)

    # Return spectrum, wavelength coordinates, and pixel size in km/s
    return specNew, LamNew, velscale



def log_rebin_pix(lamRange, n, oversample=None, velscale=None, log10=False, newRange=None):
    """
    Determine the number of new pixels and their coordinate step when
    rebinning a spectrum in geometrically stepped bins.  The input
    spectrum must be sampled linearly in wavelength.  This is primarily
    a support routine for :func:`log_rebin`.

    Although breaking this out from the main :func:`log_rebin` function
    leads to a few repeat calculations in that function, the use of this
    function is in determine a common wavelength range for a large
    number of spectra before resampling the spectra themselves.  See
    :class:`mangadap.TplLibrary` .

    Args:
        lamRange (numpy.ndarray): two elements vector containing the
            central wavelength of the first and last pixels in the
            spectrum, which is assumed to have constant wavelength
            scale! E.g. from the values in the standard FITS keywords:
            LAMRANGE = CRVAL1 + [0,CDELT1*(NAXIS1-1)].  It must be
            LAMRANGE[0] < LAMRANGE[1].
        n (int): Number of pixels in the original spectrum.
        oversample (int): (Optional) Oversampling can be done, not to
            loose spectral resolution, especally for extended wavelength
            ranges and to avoid aliasing.  Default is to provide the
            same number of output pixels as input.
        velscale (float): Velocity scale in km/s per pixels. If this
            variable is not defined, then it will contain in output the
            velocity scale.  If this variable is defined by the user it
            will be used to set the output number of pixels and
            wavelength scale.
        log10 (bool): (Optional) Flag that the spectrum should be binned
            in units of base-10 log wavelength, instead of natural log
        newRange (numpy.ndarray): (Optional) Force the spectrum to be
            sampled to a new spectral range (lamRange is the *existing*
            spectral range).

    Returns:
        float, int : Returns (1) the linear wavelength step of each
            pixel in the input spectrum, (2) the number of pixels for
            the rebinned spectrum, (3) the log-linear wavelength step
            for each pixel in the new spectrum, (4) the velocity step
            for each pixel in the new spectrum.

    Raises:
        ValueError: Raised if the input wavelength range (*lamRange* or
            *newRange*) does not have two elements or is not sorted.
        
    """
    lamRange = numpy.asarray(lamRange)
    if len(lamRange) != 2:
        raise ValueError('lamRange must contain two elements')
    if lamRange[0] >= lamRange[1]:
        raise ValueError('It must be lamRange[0] < lamRange[1]')

    # Size of output spectrum
    m = int(n) if oversample is None else int(n*oversample)

    # Get the sampling of the existing spectrum
    dLam = numpy.diff(lamRange)/(n - 1.)        # Assume constant dLam
    lim = lamRange/dLam + [-0.5, 0.5]           # All in units of dLam

    # Get the sampling for the new spectrum, if requested to be
    # different
    if newRange is not None:
        newRange = numpy.asarray(newRange)
        if len(newRange) != 2:
            raise ValueError('newRange must contain two elements')
        if newRange[0] >= newRange[1]:
            raise ValueError('It must be newRange[0] < newRange[1]')
        lim = newRange/dLam + [-0.5, 0.5]       # Set limits to a new wavelength range

        # Adjust the length
        nn = int((lamRange[1]-lamRange[0]-newRange[1]+newRange[0])/dLam)
        m = m-nn if oversample is None else m-nn*oversample

    # Set the limits to the (base-10 or natural) log of the wavelength
    logLim = numpy.log(lim) if not log10 else numpy.log10(lim)

    c = astropy.constants.c.to('km/s').value    # Speed of light in km/s (use astropy definition)

    # Set the velocity scale, if velscale not provided; otherwise force
    # the sampling based in the input velscale
    if velscale is None:                        # Velocity scale is not set by user
        velscale = numpy.diff(logLim)[0]/m*c    # Only for output
        if log10:
            velscale *= numpy.log(10.)          # Adjust to log base-10

    logscale = velscale/c                       # dlambda/lambda = dln(lambda)
    if log10:
        logscale /= numpy.log(10.)              # Convert to dlog10(lambda)
    m = int(numpy.diff(logLim)/logscale)        # Number of output pixels

    return dLam, m, logscale, velscale






