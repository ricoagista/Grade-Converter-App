
# 🎓 Program Konversi Nilai ke Huruf (Java)

Program ini dibuat untuk mengonversi **nilai angka (0–100)** menjadi **nilai huruf (A, B, C, D, E)** serta memberikan **keterangan Lulus atau Tidak Lulus**.  
Program ditulis menggunakan **Java** dan menerapkan konsep **Object-Oriented Programming (OOP)** dengan pembagian kelas yang rapi dan mudah dipahami.

---

## 🧩 Struktur Folder

```
src/
 ├── Main.java
 └── converter/
      ├── GradeConverter.java
      └── Result.java
```

---

## 📄 Penjelasan Kelas dan Fungsinya

### 1️⃣ Result.java
Kelas ini digunakan untuk menyimpan hasil konversi nilai huruf dan status kelulusan.  
Terdiri dari atribut:
- `grade`: nilai huruf (A, B, C, D, E)
- `status`: keterangan lulus atau tidak lulus

Method penting:
- `getGrade()` → mengembalikan nilai huruf
- `getStatus()` → mengembalikan status kelulusan
- `toString()` → menampilkan hasil akhir dalam bentuk teks

---

### 2️⃣ GradeConverter.java
Kelas ini berisi logika konversi nilai angka ke huruf serta menentukan keterangan kelulusan.

Kriteria konversi:
| Nilai Angka | Nilai Huruf | Keterangan     |
|--------------|-------------|----------------|
| 85 – 100     | A           | Lulus          |
| 70 – 84      | B           | Lulus          |
| 55 – 69      | C           | Lulus          |
| 40 – 54      | D           | Tidak Lulus    |
| 0 – 39       | E           | Tidak Lulus    |

Method utama:
- `convert(double nilai)` → menerima input nilai angka dan menghasilkan objek `Result`.

---

### 3️⃣ Main.java
Kelas utama program. Bertugas untuk:
1. Membaca input dari pengguna (Scanner)
2. Mengonversi nilai menggunakan kelas `GradeConverter`
3. Menampilkan hasil berupa nilai huruf dan keterangan

Contoh output:
```
Masukkan nilai (0 - 100): 75
Nilai Huruf: B | Keterangan: Lulus
```

---

## 🧭 Cara Menjalankan Program di IntelliJ IDEA

1. **Buka IntelliJ IDEA** dan pilih **New Project → Java**.
2. Buat package baru bernama `converter` di dalam folder `src`.
3. Tambahkan file berikut:
    - `Main.java` di luar package
    - `GradeConverter.java` dan `Result.java` di dalam package `converter`
4. Copy kode dari masing-masing file.
5. Klik kanan pada `Main.java` → **Run 'Main'**.

---

## 📷 Contoh Tampilan di Console

```
Masukkan nilai (0 - 100): 90
Nilai Huruf: A | Keterangan: Lulus
```
ken