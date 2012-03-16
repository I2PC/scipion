#define data structures
from ctypes import c_int, c_float, c_char_p
from ctypes import Structure
class ParticlePickingStructEmx(Structure):
    _fields_ = [("coordinate_x", c_float),
                ("coordinate_y", c_float),
                ("BlockName",c_char_p)
               ]
    prefix = "_emx_particle."
    def __init__(self, coordinate_x=-1., coordinate_y=-1.,BlockName="dummyBlock"):
        super(ParticlePickingStructEmx, self).__init__( coordinate_x, coordinate_y,BlockName)

class ParticlePickingStructXmd(Structure):
    _fields_ = [("Xcoor", c_int),
                ("Ycoor", c_int),
                ("BlockName",c_char_p)
               ]
    prefix = "_"
    
    def __init__(self, Xcoor=-1, Ycoor=-1,BlockName="dumyBlock"):
        super(ParticlePickingStructXmd, self).__init__( Xcoor, Ycoor,BlockName)
            

class CtfMicrographStructEmx(Structure):
    _fields_ = [
            ("url", c_char_p),
            ("magnification",c_float),
            ("scanner_pixel_size",c_float),
            ("defocusU",c_float),
            ("defocusV",c_float),
            ("astigmatism_angle",c_float),
            ("voltage",c_float),
            ("Cs",c_float),
            ("amplitude_contrast",c_float),
            ("BlockName",c_char_p)
            ]
    prefix = "_emx_micrograph."
    def __init__(self, 
            url="",
            magnification=1.,
            scanner_pixel_size=1.,
            defocusU=0.,
            defocusV=0.,
            astigmatism_angle=0.,
            voltage=0.,
            Cs=0.,
            amplitude_contrast=0.,
            BlockName="dummyBlock"
                 ):
        super(CtfMicrographStructEmx, self).__init__( 
            url,
            magnification,
            scanner_pixel_size,
            defocusU,
            defocusV,
            astigmatism_angle,
            voltage,
            Cs,
            amplitude_contrast,
            BlockName
            )


class CtfMicrographStructXmd(Structure):
    _fields_ = [("image", c_char_p),
                ("CTF_Sampling_rate", c_float),
                ("CTF_Defocus_U",c_float),
                ("CTF_Defocus_V",c_float),
                ("CTF_Defocus_angle",c_float),
                ("CTF_Voltage",c_float),
                ("CTF_Spherical_aberration",c_float),
                ("CTF_Q0",c_float),
                ("BlockName",c_char_p)
               ]
    prefix = "_"
    
    def __init__(self, 
            image="",
            CTF_Sampling_rate=1.,
            CTF_Defocus_U=0.,
            CTF_Defocus_V=0.,
            CTF_Defocus_angle=0.,
            CTF_Voltage=0.,
            CTF_Spherical_aberration=0.,
            CTF_Q0=0.,
            BlockName="dumyBlock"
            ):
        
        super(CtfMicrographStructXmd, self).__init__(
            image,
            CTF_Sampling_rate,
            CTF_Defocus_U,
            CTF_Defocus_V,
            CTF_Defocus_angle,
            CTF_Voltage,
            CTF_Spherical_aberration,
            CTF_Q0,
            BlockName
            )
###############Aligment
class ParticleAlignmentEMX(Structure):
    _fields_ = [
            ("url", c_char_p),
            ("magnification",c_float),
            ("scanner_pixel_size",c_float),
            ("defocusU",c_float),
            ("defocusV",c_float),
            ("astigmatism_angle",c_float),
            ("voltage",c_float),
            ("Cs",c_float),
            ("amplitude_contrast",c_float),
            ("BlockName",c_char_p)
            ]
    prefix = "_emx_micrograph."
    def __init__(self, 
            url="",
            magnification=1.,
            scanner_pixel_size=1.,
            defocusU=0.,
            defocusV=0.,
            astigmatism_angle=0.,
            voltage=0.,
            Cs=0.,
            amplitude_contrast=0.,
            BlockName="dummyBlock"
                 ):
        super(ParticleAlignmentEMX, self).__init__( 
            url,
            magnification,
            scanner_pixel_size,
            defocusU,
            defocusV,
            astigmatism_angle,
            voltage,
            Cs,
            amplitude_contrast,
            BlockName
            )


class ParticleAlignmentXmd(Structure):
    _fields_ = [("image", c_char_p),
                ("CTF_Sampling_rate", c_float),
                ("CTF_Defocus_U",c_float),
                ("CTF_Defocus_V",c_float),
                ("CTF_Defocus_angle",c_float),
                ("CTF_Voltage",c_float),
                ("CTF_Spherical_aberration",c_float),
                ("CTF_Q0",c_float),
                ("BlockName",c_char_p)
               ]
    prefix = "_"
    
    def __init__(self, 
            image="",
            CTF_Sampling_rate=1.,
            CTF_Defocus_U=0.,
            CTF_Defocus_V=0.,
            CTF_Defocus_angle=0.,
            CTF_Voltage=0.,
            CTF_Spherical_aberration=0.,
            CTF_Q0=0.,
            BlockName="dumyBlock"
            ):
        
        super(ParticleAlignmentXmd, self).__init__(
            image,
            CTF_Sampling_rate,
            CTF_Defocus_U,
            CTF_Defocus_V,
            CTF_Defocus_angle,
            CTF_Voltage,
            CTF_Spherical_aberration,
            CTF_Q0,
            BlockName
            )
##################aligment
class BlockNamesEMX(Structure):
    _fields_ = [("BlockName", c_char_p),
                ("size", c_int)]
    def __init__(self, BlockName, size=1):
        super(BlockNamesEMX, self).__init__( BlockName, size)

class BlockNamesXMD(Structure):
    _fields_ = [("BlockName", c_char_p),
                ("size", c_int)]
    def __init__(self, BlockName, size=1):
        super(BlockNamesXMD, self).__init__( BlockName, size)

