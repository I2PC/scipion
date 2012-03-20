#define data structures
from ctypes import c_int, c_float, c_char_p
from ctypes import Structure
prefix_micrograph = "_emx_micrograph."
prefix_particle   = "_emx_particle."
prefix_class      = "_emx_class."

class ParticlePickingStructEmx(Structure):
    _fields_ = [("coordinate_x", c_float),
                ("coordinate_y", c_float),
                ("url",c_char_p)
               ]
    prefix = prefix_particle
    def __init__(self, coordinate_x=-1., coordinate_y=-1.,url="dummyMic"):
        super(ParticlePickingStructEmx, self).__init__( coordinate_x, coordinate_y,url)

class ParticlePickingStructXmd(Structure):
    _fields_ = [("Xcoor", c_int),
                ("Ycoor", c_int),
                ("MicName",c_char_p)
               ]
    prefix = "_"
    
    def __init__(self, Xcoor=-1, Ycoor=-1,MicName="dummyMic"):
        super(ParticlePickingStructXmd, self).__init__( Xcoor, Ycoor,MicName)
            

class ParticleClassStructEmx(Structure):
    _fields_ = [("id", c_char_p),
                ("url",c_char_p)
               ]
    prefix = prefix_class
    def __init__(self, id="dummy", url="dummyMic"):
        super(ParticleClassStructEmx, self).__init__( id, url)

class ParticleClassStructXmd(Structure):
    _fields_ = [("image", c_char_p),
                ("ref", c_int)
               ]
    prefix = "_"
    
    def __init__(self, image="dummy.mrc", ref=-1):
        super(ParticleClassStructXmd, self).__init__( image, ref)
            

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
            ("amplitude_contrast",c_float)
            ]
    prefix = prefix_micrograph
    def __init__(self, 
            url="",
            magnification=1.,
            scanner_pixel_size=1.,
            defocusU=0.,
            defocusV=0.,
            astigmatism_angle=0.,
            voltage=0.,
            Cs=0.,
            amplitude_contrast=0.
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
            amplitude_contrast
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

class ParticleAlignmentStructEmx(Structure):
    _fields_ = [
                ("url",c_char_p),
                ("transformation_matrix_1_1",c_float),
                ("transformation_matrix_1_2",c_float),
                ("transformation_matrix_1_3",c_float),
                ("transformation_matrix_offset_x",c_float),
                ("transformation_matrix_2_1",c_float),
                ("transformation_matrix_2_2",c_float),
                ("transformation_matrix_2_3",c_float),
                ("transformation_matrix_offset_y",c_float),
                ("transformation_matrix_3_1",c_float),
                ("transformation_matrix_3_2",c_float),
                ("transformation_matrix_3_3",c_float),
                ("transformation_matrix_offset_z",c_float),
                ("enable",c_int),
                ("FOM",c_float)
                ]
    prefix = prefix_particle
    def __init__(self, 
            url="",
            transformation_matrix_1_1=1.,
            transformation_matrix_1_2=0.,
            transformation_matrix_1_3=0.,
            transformation_matrix_offset_x=0.,
            transformation_matrix_2_1=0.,
            transformation_matrix_2_2=1.,
            transformation_matrix_2_3=0.,
            transformation_matrix_offset_y=0.,
            transformation_matrix_3_1=0.,
            transformation_matrix_3_2=0.,
            transformation_matrix_3_3=1.,
            transformation_matrix_offset_z=0.,
            enable=1,
            FOM=1.
                 ):
        super(ParticleAlignmentStructEmx, self).__init__( 
                url,
                transformation_matrix_1_1,
                transformation_matrix_1_2,
                transformation_matrix_1_3,
                transformation_matrix_offset_x,
                transformation_matrix_2_1,
                transformation_matrix_2_2,
                transformation_matrix_2_3,
                transformation_matrix_offset_y,
                transformation_matrix_3_1,
                transformation_matrix_3_2,
                transformation_matrix_3_3,
                transformation_matrix_offset_z,
                enable,
                FOM
            )

class ParticleAlignmentStructXmd(Structure):
    _fields_ = [
                ("image",c_char_p),
                ("angleRot",c_float),
                ("angleTilt",c_float),
                ("anglePsi",c_float),
                ("shiftX",c_float),
                ("shiftY",c_float),
                ("shiftZ",c_float),
                ("flip",c_int),
                ("scale",c_float),
                ("enable",c_int),
                ("fom",c_float)
               ]
    prefix = "_"
    
    def __init__(self, 
            image="",
            angleRot=0.,
            angleTilt=0.,
            anglePsi=0.,
            shiftX=0.,
            shiftY=0.,
            shiftZ=0.,
            flip=0,
            scale=1.,
            enable=1,
            fom=1.
            ):
        
        super(ParticleAlignmentStructXmd, self).__init__(
            image,
            angleRot,
            angleTilt,
            anglePsi,
            shiftX,
            shiftY,
            shiftZ,
            flip,
            scale,
            enable,
            fom
            )
##################aligment_end
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

