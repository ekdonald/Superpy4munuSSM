/**
 * @file tHTables.cpp
 *
 * @copyright Copyright 2019 by the authors.
 * HiggsBounds is released under the GPLv3+.
 *
 * This file is part of HiggsBounds.
 * The Coefficients were extracted from a LO MG5 calculation in the Higgs
 * characterization model. The tWH cxn is the result of the same LO calculation,
 * with the exception of the 125GeV datapoint, which is the (only available) NLO
 * YR4 value.
 */
#include "utilities/ArithmeticArray.hpp"
#include <array>
#include <vector>

#pragma GCC diagnostic ignored "-Wmissing-braces"

namespace Higgs {
namespace {
const auto tHMassGrid = std::array{std::vector<double>{
    25.,  50.,  75.,  100., 125., 150., 175., 200., 225., 250.,
    275., 300., 325., 350., 375., 400., 425., 450., 475., 500.,
    525., 550., 575., 600., 625., 650., 675., 700., 725., 750.,
    775., 800., 825., 850., 875., 900., 925., 950., 975., 1000.}};

const auto tHCoeffs = std::vector<utilities::ArithmeticArray<double, 4>>{
    {2.14491978, 0.19858888, 0.58785621, -1.732776},
    {2.87093159, 0.43779872, 1.34377861, -3.2147102},
    {3.3530206, 0.71557965, 2.46302251, -4.81604311},
    {3.54288816, 0.92670093, 3.36545807, -5.90834623},
    {3.18448956, 0.98054168, 3.77087255, -5.95536211},
    {2.77184902, 0.97736048, 3.91164699, -5.68349601},
    {2.32077225, 0.91444388, 4.0051168, -5.32588905},
    {1.94213324, 0.84507999, 3.69985456, -4.64198781},
    {1.64178289, 0.747187, 3.47055888, -4.11234177},
    {1.40153061, 0.67511776, 3.26625882, -3.66778943},
    {1.19747956, 0.61134647, 3.08471548, -3.28219504},
    {1.02563959, 0.55509935, 2.92380095, -2.94944054},
    {0.88224613, 0.50564194, 2.7814983, -2.66374443},
    {0.76376015, 0.46227918, 2.65590159, -2.41966174},
    {0.66686818, 0.42435539, 2.54521588, -2.21208406},
    {0.5884823, 0.39125424, 2.44775721, -2.03623951},
    {0.52574017, 0.36239878, 2.36195261, -1.88769278},
    {0.47600499, 0.33725142, 2.2863401, -1.76234508},
    {0.4368655, 0.31531397, 2.21956869, -1.65643419},
    {0.40613603, 0.29612759, 2.16039838, -1.56653441},
    {0.38185645, 0.27927281, 2.10770017, -1.48955662},
    {0.36229219, 0.26436953, 2.06045602, -1.42274821},
    {0.34593422, 0.25107704, 2.01775892, -1.36369314},
    {0.3314991, 0.23909399, 1.9788128, -1.3103119},
    {0.31792891, 0.22815839, 1.94293264, -1.26086155},
    {0.30439132, 0.21804763, 1.90954435, -1.21393567},
    {0.29027953, 0.20857849, 1.87818487, -1.1684644},
    {0.27521232, 0.19960709, 1.8485021, -1.12371442},
    {0.259034, 0.19102895, 1.82025497, -1.07928897},
    {0.24181447, 0.18277893, 1.79331335, -1.03512782},
    {0.22384915, 0.17483128, 1.76765814, -0.99150729},
    {0.20565904, 0.16719964, 1.74338121, -0.94904024},
    {0.18799069, 0.15993698, 1.72068541, -0.9086761},
    {0.17181622, 0.15313567, 1.69988461, -0.87170082},
    {0.15833328, 0.14692745, 1.68140363, -0.83973691},
    {0.1489651, 0.14148342, 1.66577833, -0.81474342},
    {0.14536045, 0.13701405, 1.6536555, -0.79901595},
    {0.14939367, 0.1337692, 1.64579297, -0.79518665},
    {0.16316466, 0.13203809, 1.64305954, -0.80622419},
    {0.18899885, 0.1321493, 1.64643498, -0.83543383}};

const auto ttHCoeffs = std::vector<double>{
    0.04724174, 0.106674,   0.18649609, 0.29499832, 0.43565068, 0.61699045,
    0.77808344, 0.99761501, 1.18567633, 1.3255745,  1.42734578, 1.49738091,
    1.54138535, 1.56441581, 1.57091668, 1.56475656, 1.54926478, 1.52726782,
    1.50112585, 1.47276923, 1.44373496, 1.41520322, 1.38803381, 1.36280268,
    1.33983843, 1.31925875, 1.30100698, 1.28488852, 1.27060741, 1.25780276,
    1.24608527, 1.2350737,  1.22443137, 1.21390269, 1.20334959, 1.19278804,
    1.18242455, 1.17269264, 1.16428938, 1.15821179};

const auto tWHCoeffs = std::vector<utilities::ArithmeticArray<double, 4>>{
    {0.82270437, 0.15952126, 0.24481785, -0.06752223},
    {1.12010247, 0.42166413, 0.49481879, -0.61492127},
    {1.53064807, 0.79396099, 0.89904009, -1.42968815},
    {2.12559202, 1.37241787, 1.4820496, -2.60764162},
    {2.87904699, 2.1016822, 2.11508975, -3.99413674},
    {3.67074103, 2.92647134, 2.84634674, -5.51708778},
    {4.38185798, 3.86978248, 3.56169425, -6.94355223},
    {5.282422, 4.66818524, 3.98088995, -8.26331195},
    {5.65851856, 5.10390049, 4.21395372, -8.87247228},
    {5.95831858, 5.45432292, 4.32903292, -9.2873515},
    {6.119652, 5.67855122, 4.34998896, -9.46964097},
    {6.17701351, 5.80264748, 4.29779297, -9.47480648},
    {6.16489779, 5.85267379, 4.19069455, -9.35559234},
    {6.11292733, 5.85150782, 4.04439064, -9.15731797},
    {6.03123592, 5.80928962, 3.87219424, -8.90343016},
    {5.92508514, 5.7329748, 3.6852033, -8.61028843},
    {5.79973657, 5.62951899, 3.49246941, -8.29220598},
    {5.66045182, 5.50587779, 3.3011667, -7.96161852},
    {5.51249247, 5.36900684, 3.11676056, -7.62925303},
    {5.3611201, 5.22586175, 2.9431765, -7.3042966},
    {5.21159632, 5.08339814, 2.78296889, -6.99456521},
    {5.0691827, 4.94857164, 2.6374898, -6.70667251},
    {4.93809859, 4.82692545, 2.5070578, -6.44515639},
    {4.81839423, 4.71835317, 2.39112672, -6.20952095},
    {4.70907765, 4.62133599, 2.28845448, -5.99753214},
    {4.60915685, 4.5343551, 2.19727189, -5.80642875},
    {4.51763984, 4.45589167, 2.11545142, -5.63309127},
    {4.43353463, 4.3844269, 2.04067603, -5.47421066},
    {4.35584922, 4.31844198, 1.97060794, -5.32645716},
    {4.28359162, 4.25641809, 1.90305746, -5.18664908},
    {4.21576984, 4.19683642, 1.83615176, -5.0519216},
    {4.15139188, 4.13817815, 1.76850368, -4.91989556},
    {4.08946575, 4.07892447, 1.69938053, -4.78884628},
    {4.02899947, 4.01755658, 1.62887288, -4.65787235},
    {3.96900103, 3.95255565, 1.55806336, -4.52706439},
    {3.90847844, 3.88240287, 1.48919547, -4.39767391},
    {3.84643972, 3.80557943, 1.42584237, -4.27228209},
    {3.78189287, 3.72056652, 1.37307566, -4.15496853},
    {3.71384589, 3.62584533, 1.33763423, -4.05148012},
    {3.6413068, 3.51989703, 1.32809299, -3.96939979}};

const auto tWHCxnsLHC13 = std::vector<double>{
    4.0176400e-01, 1.3646430e-01, 5.9817800e-02, 3.0425500e-02, 1.7191800e-02,
    1.0786950e-02, 7.3363400e-03, 5.3866200e-03, 4.2122170e-03, 3.4745270e-03,
    2.9659580e-03, 2.5715380e-03, 2.3073270e-03, 2.1015380e-03, 1.8939700e-03,
    1.7035150e-03, 1.5605845e-03, 1.4700844e-03, 1.2993074e-03, 1.2144797e-03,
    1.1110825e-03, 1.0339067e-03, 9.5335160e-04, 8.9449130e-04, 8.3375110e-04,
    7.5831940e-04, 7.0076270e-04, 6.4437680e-04, 5.9727370e-04, 5.6029560e-04,
    5.1194521e-04, 4.8151724e-04, 4.4642053e-04, 4.1460529e-04, 3.8486033e-04,
    3.6168888e-04, 3.2963804e-04, 3.0831694e-04, 2.8875673e-04, 2.7142641e-04};

} // namespace
} // namespace Higgs