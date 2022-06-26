
Select *
From [Portfolio Project]..CovidDeaths
where continent is not null
order by 3,4

Select *
From [Portfolio Project]..CovidVaccinations	
order by 3,4

--Select Data that we are going to be using 

Select Location, date, total_cases, new_cases, total_deaths, population
From [Portfolio Project]..CovidDeaths	
order by 1,2

--Looking at Total Cases vs Total Deaths
-- Shows the likelihood of dying if you contract covid in your country ie India

Select Location, date, total_cases, total_deaths, (total_deaths/total_cases)*100 as DeathPercentage
From [Portfolio Project]..CovidDeaths	
Where Location like '%India%'
order by 1,2


--Looking at Total Cases vs Population
--Shows what percentage of population got Covid 

Select Location, date, total_cases, Population, (total_cases/Population)*100 as PercentPopulationInfected
From [Portfolio Project]..CovidDeaths	
Where Location like '%India%'
order by 1,2

--Looking at Countries with Highest Infection Rate compared to Population

Select Location, MAX(total_cases) as HighestInfectionCount, Population, MAX((total_cases/Population))*100 as PercentPopulationInfected
From [Portfolio Project]..CovidDeaths	
Where Location like '%India%'
Group by Location, Population
order by PercentPopulationInfected desc


--Showing the Countries with Highest Death Count per Population

Select Location, MAX(cast(Total_deaths as int)) as TotalDeathCount
From [Portfolio Project]..CovidDeaths	
Where Location like '%India%'
where continent is not null
Group by Location
order by TotalDeathCount desc


--LETS BREAK THINGS DOWN BY CONTINENT 

Select continent, MAX(cast(Total_deaths as int)) as TotalDeathCount
From [Portfolio Project]..CovidDeaths	
Where Location like '%India%'
where continent is not null
Group by continent
order by TotalDeathCount desc

-- Showing Continents with the highest death count per population
Select continent, MAX(cast(Total_deaths as int)) as TotalDeathCount
From [Portfolio Project]..CovidDeaths	
Where Location like '%India%'
where continent is not null
Group by continent
order by TotalDeathCount desc 

--GLOBAL NUMBERS

Select SUM(new_cases) as total_cases, SUM(cast(new_deaths as int)) as total_death
, SUM(cast(new_deaths as int))/SUM(New_cases)*100 as DeathPercentage
From [Portfolio Project]..CovidDeaths	
Where Location like '%India%'
where continent is not null
Group by date
order by 1,2

--Looking at Total Population vs Vaccinations

Select dea.continent, dea.location, dea.date, dea.population,vac.new_vaccinations
, SUM(Convert(bigint,vac.new_vaccinations)) OVER (Partition by dea.location Order by dea.location, dea.date) as Rollingpeople_vaccinated
From [Portfolio Project]..CovidDeaths dea
Join [Portfolio Project]..CovidVaccinations vac
    on dea.location = vac.location
    and dea.date = vac.date
where dea.continent is not null
--order by 2,3

--USE CTE

with PopvsVac (Continent, location, date, Population, New_vaccinations, RollingPeople_vaccinated) 
as 
(
Select dea.continent, dea.location, dea.date, dea.population,vac.new_vaccinations
, SUM(Convert(bigint,vac.new_vaccinations)) OVER (Partition by dea.location Order by dea.location, dea.date) as Rollingpeople_vaccinated
,(RollingPeopleVaccinated/Population)*100
From [Portfolio Project]..CovidDeaths dea
Join [Portfolio Project]..CovidVaccinations vac
    on dea.location = vac.location
    and dea.date = vac.date
where dea.continent is not null
order by 2,3
)
Select*, (RollingPeople_vaccinated/Population)*100 as Vaccinatedpercentage
From PopvsVac


--TEMP TABLE
DROP TABLE if exists #PercentPopulationVaccinated
Create Table #PercentPopulationVaccinated
(
Continent nvarchar(255),
location nvarchar(255),
Date datetime,
Population numeric,
New_vaccinations numeric,
Rollingpeople_vaccinated numeric
)

Insert into #PercentPopulationVaccinated
Select dea.continent, dea.location, dea.date, dea.population, vac.new_vaccinations
, SUM(Convert(bigint,vac.new_vaccinations)) OVER (Partition by dea.location Order by dea.location, dea.date) as Rollingpeople_vaccinated
,(RollingPeopleVaccinated/Population)*100
From [Portfolio Project]..CovidDeaths dea
Join [Portfolio Project]..CovidVaccinations vac
   on dea.location = vac.location
   and dea.date = vac.date
where dea.continent is not null
order by 2,3

Select*, (RollingPeople_vaccinated/Population)*100 as Vaccinatedpercentage
From #PercentPopulationVaccinated


--Creating view to store data for later visualization

Create View PercentPopulationVaccinated as
Select dea.continent, dea.location, dea.date, dea.population, vac.new_vaccinations
, SUM(Convert(bigint,vac.new_vaccinations)) OVER (Partition by dea.location Order by dea.location, dea.date) as Rollingpeople_vaccinated
--,(RollingPeopleVaccinated/Population)*100
From [Portfolio Project]..CovidDeaths dea
Join [Portfolio Project]..CovidVaccinations vac
   on dea.location = vac.location
   and dea.date = vac.date
where dea.continent is not null
--order by 2,3

Select * 
From PercentPopulationVaccinated
