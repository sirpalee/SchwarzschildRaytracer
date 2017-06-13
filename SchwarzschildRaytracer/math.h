#pragma once

static constexpr double Pi = 3.1415926535897932384626433;
static constexpr double deg = Pi / 180.0;

auto sq = [](auto const& x) { return x * x; };
auto cube = [](auto const& x) { return x * x * x; };
