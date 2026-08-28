# lua with complex numbers
add_library(luacomplex
    liblua/femmcomplex.cpp
    liblua/llimits.h
    liblua/lopcodes.h
    liblua/lua.h
    liblua/luadebug.h
    liblua/lualib.h
    liblua/lapi.cpp
    liblua/lauxlib.cpp
    liblua/lbaselib.cpp
    liblua/lcode.cpp
    liblua/ldblib.cpp
    liblua/ldebug.cpp
    liblua/ldo.cpp
    liblua/lfunc.cpp
    liblua/lgc.cpp
    liblua/liolib.cpp
    liblua/llex.cpp
    liblua/lmathlib.cpp
    liblua/lmem.cpp
    liblua/lobject.cpp
    liblua/lparser.cpp
    liblua/lstate.cpp
    liblua/lstring.cpp
    liblua/lstrlib.cpp
    liblua/ltable.cpp
    liblua/ltests.cpp
    liblua/ltm.cpp
    liblua/lundump.cpp
    liblua/lvm.cpp
    liblua/lzio.cpp
    )
target_include_directories(luacomplex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/liblua> $<INSTALL_INTERFACE:include/liblua>)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/femmversion.h.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/femmversion.h"
    @ONLY
    )

add_library(femm
    femmconstants.cpp
    femmenums.cpp
    CArcSegment.cpp
    CBlockLabel.cpp
    CBoundaryProp.cpp
    CCircuit.cpp
    CCommonPoint.cpp
    CElement.cpp
    CAirGapElement.cpp
    CliTools.cpp
    CMaterialProp.cpp
    CMeshNode.cpp
    CNode.cpp
    CPointProp.cpp
    CSegment.cpp
    cspars.cpp
    cuthill.cpp
    feasolver.cpp
    FemmProblem.cpp
    FemmReader.cpp
    FemmStateBase.cpp
    femmversion.cpp
    fparse.cpp
    fullmatrix.cpp
    IntPoint.cpp
    locationTools.cpp
    LuaInstance.cpp
    MatlibReader.cpp
    PostProcessor.cpp
    spars.cpp
    stringTools.cpp
    )
target_include_directories(femm
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
target_link_libraries(femm PUBLIC luacomplex)

# Optional Eigen-based sparse direct solver.  If a copy of Eigen is
# present in cfemm/external, spars.cpp gains an LDLT direct-solve path,
# used by default (set XFEMM_DIRECT=0 in the environment to force the
# legacy PCG solver).  EIGEN_MPL2_ONLY guarantees at compile time that
# only MPL2-licensed Eigen code is included (any LGPL-licensed file
# would trigger a compile error).
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../external/eigen-3.4.0/Eigen")
    message(STATUS "Eigen found in cfemm/external - enabling XFEMM_HAVE_EIGEN direct solver")
    target_include_directories(femm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../external/eigen-3.4.0")
    target_compile_definitions(femm PRIVATE XFEMM_HAVE_EIGEN EIGEN_MPL2_ONLY)
endif()
# vi:expandtab:tabstop=4 shiftwidth=4:
