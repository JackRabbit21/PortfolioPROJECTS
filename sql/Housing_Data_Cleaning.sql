Select *
From Housing

--Standardize the sale date format

Select SaleDateConverted,Convert(Date,SaleDate)
From Housing

Update Housing
SET SaleDate = Convert(Date,SaleDate)

ALTER TABLE Housing
Add SaleDateConverted Date;

Update Housing
SET SaleDateConverted = Convert(Date,SaleDate)

--Populate Property Address data

Select *
From Housing
--Where PropertyAddress is null
order by ParcelID


Select a.ParcelID,a.PropertyAddress,b.ParcelID,b.PropertyAddress,ISNULL(a.propertyAddress,b.PropertyAddress)
From [Portfolio Project]..Housing a
JOIN [Portfolio Project]..Housing b
    on a.ParcelID = b.ParcelID
	AND a.[UniqueID ]<> b.[UniqueID ]
Where a.PropertyAddress is null

Update a
SET PropertyAddress = ISNULL(a.propertyaddress,b.PropertyAddress)
From [Portfolio Project]..Housing a
JOIN [Portfolio Project]..Housing b
    on a.ParcelID = b.ParcelID
	AND a.[UniqueID ]<> b.[UniqueID ]

-- Breaking out Address into Individual Columns (Address, City, State)

Select PropertyAddress
From [Portfolio Project]..Housing
--Where PropertyAddress is null
--order by ParcelID

Select
SUBSTRING(PropertyAddress, 1,CHARINDEX(',' , PropertyAddress)-1) as Address
,SUBSTRING(PropertyAddress, CHARINDEX(',' , PropertyAddress) +1, LEN(PropertyAddress)) as Address
From [Portfolio Project]..Housing

ALTER TABLE Housing
Add PropertySplitAddress nvarchar(255);

Update Housing
SET PropertySplitAddress = SUBSTRING(PropertyAddress, 1,CHARINDEX(',' , PropertyAddress)-1)

ALTER TABLE Housing
Add PropertySplitCity nvarchar(255);

Update Housing
SET PropertySplitCity = SUBSTRING(PropertyAddress, CHARINDEX(',' , PropertyAddress) +1, LEN(PropertyAddress))

Select *
From [Portfolio Project]..Housing

Select OwnerAddress
From [Portfolio Project]..Housing

Select
PARSENAME(REPLACE(OwnerAddress,',','.'),3)
,PARSENAME(REPLACE(OwnerAddress,',','.'),2)
,PARSENAME(REPLACE(OwnerAddress,',','.'),1)
From [Portfolio Project]..Housing

ALTER TABLE Housing
Add OwnerSplitAddress nvarchar(255);

Update Housing
SET OwnerSplitAddress = PARSENAME(REPLACE(OwnerAddress,',','.'),3)

ALTER TABLE Housing
Add OwnerSplitCity nvarchar(255);

Update Housing
SET OwnerSplitCity = PARSENAME(REPLACE(OwnerAddress,',','.'),2)

ALTER TABLE Housing
Add OwnerSplitState nvarchar(255);

Update Housing
SET OwnerSplitState = PARSENAME(REPLACE(OwnerAddress,',','.'),1)

Select *
From Housing

--Change Y and N to Yes and No in "Sold as Vacant" field

Select Distinct(SoldASVacant),Count(SoldAsVacant)
From [Portfolio Project]..Housing
Group by SoldAsVacant
order by 2

Select SoldAsVacant
,CASE When SoldAsVacant = 'Y' THEN 'Yes'
      When SoldAsVacant = 'N' THEN 'No'
	  ELSE SoldAsVacant 
	  END
From [Portfolio Project]..Housing

Update [Portfolio Project]..Housing
SET SoldAsVacant = CASE When SoldAsVacant = 'Y' THEN 'Yes'
      When SoldAsVacant = 'N' THEN 'No'
	  ELSE SoldAsVacant 
	  END 

--Remove Duplicates
WITH RowNumCTE AS(
Select *,
   ROW_NUMBER() OVER(
   PARTITION BY ParcelID,
                PropertyAddress,
				SalePrice,
				SaleDate,
				LegalReference
				ORDER BY
				   UniqueID
				   ) row_num
		
From [Portfolio Project]..Housing
--order by ParcelID
)
Select *
From RowNumCTE
Where row_num > 1
--Order by PropertyAddress

--Delete Unused Columns

Select *
From [Portfolio Project]..Housing

AlTER TABLE[Portfolio Project]..Housing
Drop COLUMN OwnerAddress, TaxDistrict, PropertyAddress 

AlTER TABLE[Portfolio Project]..Housing
Drop COLUMN SaleDate